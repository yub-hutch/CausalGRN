#' Perform conditional independence tests on graph edges based on wild-type data
#'
#' Applies conditional independence (CI) tests to each edge in the graph `G`. Edges are evaluated in parallel,
#' and those found conditionally independent given a subset of variables are removed.
#' The function updates the graph and associated statistics accordingly.
#'
#' @param G Adjacency matrix.
#' @param order Integer specifying the size of the conditioning set.
#' @param max_order Maximum allowed size for the conditioning set.
#' @param hub_index Vector of node indices considered as hubs.
#' @param C Correlation matrix of the variables.
#' @param n Sample size used in the CI tests.
#' @param alpha Significance level for the CI tests.
#' @param pMax Matrix to store the maximum p-values for each edge.
#' @param chisqMin Matrix to store the minimum chi-square statistics for each edge.
#' @param sepSet List of separation sets for each pair of nodes.
#' @param ncores Number of cores to use for parallel processing.
#'
#' @return A list containing:
#' \describe{
#'   \item{G}{Updated adjacency matrix after removing conditionally independent edges.}
#'   \item{pMax}{Updated matrix of maximum p-values for each edge.}
#'   \item{chisqMin}{Updated matrix of minimum chi-square statistics for each edge.}
#'   \item{sepSet}{Updated list of separation sets for each pair of nodes.}
#'   \item{done}{Logical indicating whether all edges reached the specified order.}
#' }
#' @export
perform_ci_test_from_wildtype <- function(G, order, max_order, hub_index, C, n, alpha, pMax, chisqMin, sepSet, ncores) {
  edge_index <- get_edge_index(G)
  message('Performing CI test ...')
  res <- pbmcapply::pbmcmapply(function(i, j) {
    if (order == 0) {
      A <- integer(0)
    } else {
      A <- setdiff(which(G[, i] | G[, j]), c(i, j))
    }
    if (order > max_order) {
      A <- intersect(A, hub_index)
    }
    if (length(A) < order) {
      return(list(order_reached = TRUE))
    }
    pmax <- -1
    chisqmin <- Inf
    iter <- iterpc::iterpc(length(A), order)
    index <- iterpc::getnext(iter)
    repeat {
      S <- A[index]
      zstat <- pcalg::zStat(x = i, y = j, S = S, C = C, n = n)
      chisq <- zstat ^ 2
      pv <- pchisq(chisq, df = 1, lower.tail = FALSE)
      pmax <- max(pmax, pv)
      chisqmin <- min(chisqmin, chisq)
      if (pmax > alpha) {
        return(list(order_reached = FALSE, pmax = pmax, chisqmin = chisqmin, sepset = S))
      }
      index = iterpc::getnext(iter)
      if (is.null(index)) {
        return(list(order_reached = FALSE, pmax = pmax, chisqmin = chisqmin))
      }
    }
  }, i = edge_index[, 1], j = edge_index[, 2], SIMPLIFY = F, mc.cores = ncores)
  message('Updating graph ...')
  pb <- pbmcapply::progressBar(min = 0, max = length(res))
  for (k in seq_along(res)) {
    setTxtProgressBar(pb = pb, value = k)
    i <- edge_index[k, 1]
    j <- edge_index[k, 2]
    if (!(res[[k]]$order_reached)) {
      pMax[i, j] <- pMax[j, i] <- max(pMax[i, j], res[[k]]$pmax)
      chisqMin[i, j] <- chisqMin[j, i] <- min(chisqMin[i, j], res[[k]]$chisqmin)
      if (res[[k]]$pmax > alpha) {
        G[i, j] <- G[j, i] <- FALSE
        sepSet[[i]][[j]] <- sepSet[[j]][[i]] <- res[[k]]$sepset
      }
    }
  }
  close(pb)
  done <- all(sapply(res, function(x) x$order_reached))
  return(list(G = G, pMax = pMax, chisqMin = chisqMin, sepSet = sepSet, done = done))
}


#' Infer graph skeleton via conditional independence tests based on wild-type data
#'
#' Constructs graph skeleton by performing iterative conditional independence (CI) tests on gene expression data.
#' Edges are removed based on statistical independence, optionally considering hub genes with higher-order conditioning.
#'
#' @param Y scRNA-seq matrix (samples × genes).
#' @param alpha Significance level for CI tests.
#' @param ncores Number of CPU cores for parallel processing.
#' @param G Optional initial adjacency matrix; defaults to a fully connected graph without self-loops.
#' @param max_order Maximum conditioning set size for non-hub genes.
#' @param hub_genes Optional character vector of hub gene names.
#' @param max_order_hub Maximum conditioning set size for hub genes.
#'
#' @return A list with:
#' \describe{
#'   \item{graph}{An igraph object representing the inferred skeleton with edge attributes `pMax` and `chisqMin`.}
#'   \item{sepSet}{List of separation sets for each node pair.}
#' }
#'
#' @examples
#' # No hub genes
#' # x -> y -> z
#' set.seed(123)
#' n <- 1e3
#' x <- rnorm(n)
#' y <- rnorm(n, mean = x, sd = 1)
#' z <- rnorm(n, mean = y, sd = 1)
#' data <- cbind(x = x, y = y, z = z)
#' skel <- infer_skeleton_from_wildtype(
#'   Y = data,
#'   alpha = 0.05,
#'   ncores = 2,
#'   max_order = 1
#' )
#' print(igraph::as_data_frame(skel$graph))
#' plot(skel$graph)
#' get_sepset(skel, 'x', 'z')
#'
#' # With hub genes
#' # Hub genes PC 1-3, x_i -> PC_j -> z
#' set.seed(123)
#' n <- 1e3
#' X <- matrix(rnorm(5 * n), n, 5)
#' B <- cbind(
#'   c(0.4082, 0.4082, 0.4082, 0.4082, 0.4082),
#'   c(0.6325, 0.3162, 0.0000, -0.3162, -0.6325),
#'   c(0.0000, 0.7071, -0.7071, 0.0000, 0.0000)
#' )
#' PC <- X %*% B + matrix(rnorm(3 * n), n, 3)
#' print(cor(PC))
#' z <- PC %*% rep(1, 3) + rnorm(n)
#' data <- cbind(X, PC, z)
#' colnames(data) <- c(paste0('x', 1:5), paste0('PC', 1:3), 'z')
#' G <- matrix(TRUE, ncol(data), ncol(data), dimnames = list(colnames(data), colnames(data)))
#' diag(G) <- FALSE
#' G[paste0('PC', 1:3), paste0('PC', 1:3)] <- FALSE
#' skel <- infer_skeleton_from_wildtype(
#'   Y = data,
#'   alpha = 1e-3,
#'   ncores = 2,
#'   max_order = 1,
#'   hub_genes = paste0('PC', 1:3),
#'   max_order_hub = 3
#' )
#' plot(skel$graph, layout = igraph::layout_with_kk)
#' @export
infer_skeleton_from_wildtype <- function(Y, alpha, ncores, G = NULL, max_order = 1, hub_genes = NULL, max_order_hub = NULL) {
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
  # Check hub genes
  if (!is.null(hub_genes)) {
    stopifnot(all(hub_genes %in% genes))
    stopifnot(max_order_hub > max_order)
    hub_genes <- unique(hub_genes)
    hub_index <- match(hub_genes, genes)
  } else {
    max_order_hub <- max_order
  }
  # Initialize CI test statistics
  message('------------------------------------------------')
  message('Initializing CI test statistics ...')
  n <- nrow(Y)
  C <- parallel_cor(Y, ncores = ncores)
  pMax <- chisqMin <- matrix(NA, length(genes), length(genes), dimnames = list(genes, genes))
  pMax[G] <- -1
  chisqMin[G] <- Inf
  sepSet <- lapply(seq_along(genes), function(j) vector('list', length(genes)))
  # Perform CI test
  order <- 0
  while (any(G) && (order <= max_order_hub)) {
    message('------------------------------------------------')
    message(paste0('Order = ', order))
    res <- perform_ci_test_from_wildtype(
      G = G,
      order = order,
      max_order = max_order,
      hub_index = hub_index,
      C = C,
      n = n,
      alpha = alpha,
      pMax = pMax,
      chisqMin = chisqMin,
      sepSet = sepSet,
      ncores = ncores
    )
    G <- res$G
    pMax <- res$pMax
    chisqMin <- res$chisqMin
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
  graph <- adj2igraph(G = G, pMax = pMax, chisqMin = chisqMin)
  return(list(graph = graph, sepSet = sepSet))
}
