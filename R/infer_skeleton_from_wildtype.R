#' Perform conditional independence tests on graph edges based on wild-type data
#'
#' Applies conditional independence (CI) tests to each edge in the graph `G`. Edges are evaluated in parallel,
#' and those found conditionally independent given a subset of variables are removed.
#' The function updates the graph and associated statistics accordingly.
#'
#' @param G Adjacency matrix.
#' @param order Integer specifying the size of the conditioning set.
#' @param max_order Maximum allowed size for the conditioning set.
#' @param count scRNA-seq count matrix (cells x genes).
#' @param Y scRNA-seq normalized expression matrix (cells × genes). For example, log1p(total UMI corrected count).
#' @param max_thr Maximum threshold for conditional variable.
#' @param min_n1 Minimum number of samples satisfying Yk > selected threshold.
#' @param min_n2 Minimum number of samples satisfying Yk > selected threshold, Yi > 0, and Yj > 0.
#' @param alpha Significance level for the CI tests.
#' @param min_abspcor Minimum absolute value of partial correlation for kept edges.
#' @param pMax Matrix to store the maximum p-values for each edge.
#' @param chisqMin Matrix to store the minimum chi-square statistics for each edge.
#' @param absPcorMin Matrix to store the minimum absolute value of partial correlation for each edge.
#' @param sepSet List of separation sets for each pair of nodes.
#' @param ncores Number of cores to use for parallel processing.
#'
#' @return A list containing:
#' \describe{
#'   \item{G}{Updated adjacency matrix after removing conditionally independent edges.}
#'   \item{pMax}{Updated matrix of maximum p-values for each edge.}
#'   \item{chisqMin}{Updated matrix of minimum chi-square statistics for each edge.}
#'   \item{absPcorMin}{Updated matrix of minimum absolute value of partial correlation for each edge.}
#'   \item{sepSet}{Updated list of separation sets for each pair of nodes.}
#'   \item{done}{Logical indicating whether all edges reached the specified order.}
#' }
#' @export
perform_ci_test_from_wildtype <- function(
    G, order, max_order, count, Y, max_thr, min_n1, min_n2, alpha, min_abspcor, pMax, chisqMin, absPcorMin, sepSet, ncores
) {
  edge_index <- get_edge_index(G)
  message('Performing CI test ...')
  res <- pbmcapply::pbmcmapply(function(i, j) {
    if (order == 0) {
      A <- integer(0)
    } else {
      A <- setdiff(which(G[, i] | G[, j]), c(i, j))
    }
    if (length(A) < order) {
      return(list(order_reached = TRUE))
    }
    pmax <- -1
    chisqmin <- Inf
    abspcormin <- Inf
    iter <- iterpc::iterpc(length(A), order)
    index <- iterpc::getnext(iter)
    repeat {
      S <- A[index]
      if (order == 0) {
        res_test <- stats::cor.test(Y[, i], Y[, j])
      }
      if (order == 1) {
        res_test <- calc_pcor(i = i, j = j, k = S, count = count, Y = Y, max_thr = max_thr, min_n1 = min_n1, min_n2 = min_n2)
      }
      chisq <- res_test$statistic ^ 2
      pv <- res_test$p.value
      abspcor <- abs(res_test$estimate)
      pmax <- max(pmax, pv)
      chisqmin <- min(chisqmin, chisq)
      abspcormin <- min(abspcormin, abspcor)
      if (pmax > alpha || abspcormin < min_abspcor) {
        return(list(order_reached = FALSE, pmax = pmax, chisqmin = chisqmin, abspcormin = abspcormin, sepset = S))
      }
      index = iterpc::getnext(iter)
      if (is.null(index)) {
        return(list(order_reached = FALSE, pmax = pmax, chisqmin = chisqmin, abspcormin = abspcormin))
      }
    }
  }, i = edge_index[, 1], j = edge_index[, 2], SIMPLIFY = F, mc.cores = ncores)
  message('Updating graph ...')
  needs_update <- !vapply(res, `[[`, logical(1), "order_reached")
  if (any(needs_update)) {
    i <- edge_index[needs_update, 1]
    j <- edge_index[needs_update, 2]
    p <- vapply(res[needs_update], `[[`, numeric(1), "pmax")
    chi <- vapply(res[needs_update], `[[`, numeric(1), "chisqmin")
    abspcor <- vapply(res[needs_update],  `[[`, numeric(1), "abspcormin")
    idx <- cbind(i, j)
    mirror_idx <- idx[, 2:1, drop = FALSE]
    pMax[idx] <- pMax[mirror_idx] <- pmax(pMax[idx], p)
    chisqMin[idx] <- chisqMin[mirror_idx] <- pmin(chisqMin[idx], chi)
    absPcorMin[idx] <- absPcorMin[mirror_idx] <- pmin(absPcorMin[idx], abspcor)
    to_remove <- (p > alpha) | (abspcor < min_abspcor)
    if (any(to_remove)) {
      rem <- idx[to_remove, , drop = FALSE]
      mirror_rem  <- rem[, 2:1, drop = FALSE]
      G[rem] <- G[mirror_rem] <- FALSE
      if (!is.null(sepSet)) {
        pos <- which(needs_update)[to_remove]
        pb <- pbmcapply::progressBar(min = 0, max = length(pos))
        for (k in seq_along(pos)) {
          setTxtProgressBar(pb, value = k)
          a   <- rem[k, 1]
          b   <- rem[k, 2]
          sep <- res[[pos[k]]]$sepset
          sepSet[[a]][[b]] <- sepSet[[b]][[a]] <- sep
        }
        close(pb)
      }
    }
  }
  done <- !any(needs_update)
  return(list(G = G, pMax = pMax, chisqMin = chisqMin, absPcorMin = absPcorMin, sepSet = sepSet, done = done))
}


#' Infer graph skeleton via conditional independence tests based on wild-type data
#'
#' Constructs graph skeleton by performing iterative conditional independence (CI) tests on gene expression data.
#' Edges are removed based on statistical independence.
#'
#' @param count scRNA-seq count matrix (cells x genes).
#' @param Y scRNA-seq normalized expression matrix (cells × genes). For example, log1p(total UMI corrected count).
#' @param alpha Significance level for CI tests.
#' @param min_abspcor Minimum absolute value of partial correlation for kept edges.
#' @param ncores Number of CPU cores for parallel processing.
#' @param G Optional initial adjacency matrix (defaults to a fully connected graph without self-loops).
#' @param max_order Maximum conditioning set size (0 or 1, default is 1).
#' @param max_thr Maximum threshold for conditional variable (default is 10).
#' @param min_n1 Minimum number of samples satisfying Yk > selected threshold (default is 2000).
#' @param min_n2 Minimum number of samples satisfying Yk > selected threshold, Yi > 0, and Yj > 0 (default is 400).
#' @param max_nchildren Maximum number of children a node can have (default is Inf).
#' @param max_nparent Maximum number of parents a node can have (default is Inf).
#' @param sepset Return separation set or not (default is \code{TRUE}).
#'
#' @return A list with:
#' \describe{
#'   \item{graph}{An igraph object representing the inferred skeleton with edge attributes `pMax`, `chisqMin`, and `absPcorMin`.}
#'   \item{sepSet}{List of separation sets for each node pair.}
#' }
#'
#' @examples
#' # Gene 1 -> Gene 2 -> Gene 3
#' set.seed(123)
#' n <- 1e5
#' u1 <- rnorm(n)
#' u2 <- rnorm(n, u1)
#' u3 <- rnorm(n, u2)
#' u <- cbind(g1 = u1, g2 = u2, g3 = u3)
#' count <- apply(u, 2, function(x) rpois(n, lambda = exp(x)))
#' Y <- log1p(count)
#' skel <- infer_skeleton_from_wildtype(
#'   count = count,
#'   Y = Y,
#'   alpha = 1e-3,
#'   min_abspcor = 0.05,
#'   ncores = 2,
#'   max_order = 1
#' )
#' print(igraph::as_data_frame(skel$graph))
#' plot(skel$graph)
#' get_sepset(skel, 'g1', 'g3')
#'
#' @export
infer_skeleton_from_wildtype <- function(
    count, Y, alpha, min_abspcor, ncores, G = NULL, max_order = 1, max_thr = 10, min_n1 = 2000, min_n2 = 400,
    max_nchildren = Inf, max_nparent = Inf, sepset = TRUE
) {
  stopifnot(identical(dimnames(count), dimnames(Y)))
  # Check initial adjacency matrix
  genes = colnames(Y)
  if (is.null(G)) {
    G <- matrix(TRUE, length(genes), length(genes), dimnames = list(genes, genes))
    diag(G) <- FALSE
  } else {
    stopifnot(setequal(rownames(G), genes) && setequal(colnames(G), genes))
    stopifnot(!any(diag(G)))
    G <- G[genes, genes]
  }
  last_nedge <- nedge <- sum(G) / 2
  message(paste0("Number of edges in initial graph: ", nedge))
  # Initialize CI test statistics
  pMax <- chisqMin <- absPcorMin <- matrix(NA, length(genes), length(genes), dimnames = list(genes, genes))
  pMax[G] <- -1
  chisqMin[G] <- Inf
  absPcorMin[G] <- Inf
  if (sepset) {
    sepSet <- lapply(seq_along(genes), function(j) vector('list', length(genes)))
  } else {
    sepSet <- NULL
  }
  # Perform CI test
  order <- 0
  while (any(G) && (order <= max_order)) {
    message('------------------------------------------------')
    message(paste0('Order = ', order))
    res <- perform_ci_test_from_wildtype(
      G = G,
      order = order,
      max_order = max_order,
      count = count,
      Y = Y,
      max_thr = max_thr,
      min_n1 = min_n1,
      min_n2 = min_n2,
      alpha = alpha,
      min_abspcor = min_abspcor,
      pMax = pMax,
      chisqMin = chisqMin,
      absPcorMin = absPcorMin,
      sepSet = sepSet,
      ncores = ncores
    )
    G <- res$G
    pMax <- res$pMax
    chisqMin <- res$chisqMin
    absPcorMin <- res$absPcorMin
    sepSet <- res$sepSet
    done <- res$done
    nedge <- sum(G) / 2
    message(paste0('Number of edges removed: ', last_nedge - nedge))
    last_nedge <- nedge
    if (res$done) {
      message('CI tests are finished before reaching specified max order.')
      break
    }
    rm(res)
    gc()
    order <- order + 1
  }
  graph <- adj2igraph(
    G = G, pMax = pMax, chisqMin = chisqMin, absPcorMin = absPcorMin,
    max_nchildren = max_nchildren, max_nparent = max_nparent
  )
  return(list(graph = graph, sepSet = sepSet))
}
