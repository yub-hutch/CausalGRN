# Run one order of conditional independence tests and update skeleton statistics.
.perform_ci_test <- function(
    G, order, count, Y, max_thr, min_n1, min_n2,
    alpha, min_abspcor, pMax, chisqMin, absPcorMin,
    Threshold, sampleSize, sepSet, ncores
) {
  edge_index <- get_edge_index(G)
  message('Performing CI test ...')
  if (order == 0) {
    # Order-0 tests are ordinary Pearson marginal tests, so all correlations
    # can be computed once and converted to the same t statistic as cor.test().
    cors <- WGCNA::cor(Y, use = "all.obs", method = "pearson", nThreads = ncores)
    idx <- edge_index
    mirror_idx <- idx[, 2:1, drop = FALSE]
    rho <- cors[idx]
    df <- nrow(Y) - 2L
    statistic <- rho * sqrt(df / (1 - rho ^ 2))
    p <- 2 * stats::pt(abs(statistic), df = df, lower.tail = FALSE)
    chi <- statistic ^ 2
    abspcor <- abs(rho)

    message('Updating graph ...')
    pMax[idx] <- pMax[mirror_idx] <- pmax(pMax[idx], p)
    chisqMin[idx] <- chisqMin[mirror_idx] <- pmin(chisqMin[idx], chi)
    to_update <- (abspcor < absPcorMin[idx])
    absPcorMin[idx][to_update] <- absPcorMin[mirror_idx][to_update] <- abspcor[to_update]
    # Threshold and sample size only exist for order-1 partial-correlation tests.
    Threshold[idx][to_update] <- Threshold[mirror_idx][to_update] <- NA_real_
    sampleSize[idx][to_update] <- sampleSize[mirror_idx][to_update] <- NA_real_
    to_remove <- (p > alpha) | (abspcor < min_abspcor)
    if (any(to_remove)) {
      rem <- idx[to_remove, , drop = FALSE]
      mirror_rem <- rem[, 2:1, drop = FALSE]
      G[rem] <- G[mirror_rem] <- FALSE
      if (!is.null(sepSet)) {
        pb <- pbmcapply::progressBar(min = 0, max = nrow(rem))
        for (k in seq_len(nrow(rem))) {
          setTxtProgressBar(pb, value = k)
          a <- rem[k, 1]
          b <- rem[k, 2]
          sepSet[[a]][[b]] <- sepSet[[b]][[a]] <- integer(0)
        }
        close(pb)
      }
    }
    # Order 0 never finishes early due to missing conditioning sets.
    done <- FALSE
    return(list(
      G = G, pMax = pMax, chisqMin = chisqMin,
      absPcorMin = absPcorMin, Threshold = Threshold,
      sampleSize = sampleSize, sepSet = sepSet, done = done
    ))
  }

  res <- .parallel_lapply(
    seq_len(nrow(edge_index)),
    function(pos) {
      i <- edge_index[pos, 1]
      j <- edge_index[pos, 2]
      A <- setdiff(which(G[, i] | G[, j]), c(i, j))
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
        res_test <- calc_pcor(
          i = i, j = j, k = S, count = count, Y = Y,
          max_thr = max_thr, min_n1 = min_n1, min_n2 = min_n2
        )
        curr_threshold <- res_test$threshold
        curr_sample_size <- res_test$n
        chisq <- res_test$statistic ^ 2
        pv <- res_test$p.value
        abspcor <- abs(res_test$estimate)
        pmax <- max(pmax, pv)
        chisqmin <- min(chisqmin, chisq)
        if (abspcor < abspcormin) {
          abspcormin <- abspcor
          threshold <- curr_threshold
          sample_size <- curr_sample_size
        }
        if (pmax > alpha || abspcormin < min_abspcor) {
          return(list(
            order_reached = FALSE, pmax = pmax, chisqmin = chisqmin, abspcormin = abspcormin,
            threshold = threshold, sample_size = sample_size, sepset = S
          ))
        }
        index <- iterpc::getnext(iter)
        if (is.null(index)) {
          return(list(
            order_reached = FALSE, pmax = pmax, chisqmin = chisqmin, abspcormin = abspcormin,
            threshold = threshold, sample_size = sample_size
          ))
        }
      }
    },
    ncores = ncores,
    preschedule = FALSE,
    export = c(
      "edge_index", "G", "order", "count", "Y", "max_thr",
      "min_n1", "min_n2", "alpha", "min_abspcor", "calc_pcor"
    )
  )
  message('Updating graph ...')
  needs_update <- !vapply(res, `[[`, logical(1), "order_reached")
  if (any(needs_update)) {
    i <- edge_index[needs_update, 1]
    j <- edge_index[needs_update, 2]
    p <- vapply(res[needs_update], `[[`, numeric(1), "pmax")
    chi <- vapply(res[needs_update], `[[`, numeric(1), "chisqmin")
    abspcor <- vapply(res[needs_update], `[[`, numeric(1), "abspcormin")
    thr <- vapply(res[needs_update], `[[`, numeric(1), "threshold")
    ss <- vapply(res[needs_update], `[[`, numeric(1), "sample_size")
    idx <- cbind(i, j)
    mirror_idx <- idx[, 2:1, drop = FALSE]
    pMax[idx] <- pMax[mirror_idx] <- pmax(pMax[idx], p)
    chisqMin[idx] <- chisqMin[mirror_idx] <- pmin(chisqMin[idx], chi)
    to_update <- (abspcor < absPcorMin[idx])
    absPcorMin[idx][to_update] <- absPcorMin[mirror_idx][to_update] <- abspcor[to_update]
    Threshold[idx][to_update] <- Threshold[mirror_idx][to_update] <- thr[to_update]
    sampleSize[idx][to_update] <- sampleSize[mirror_idx][to_update] <- ss[to_update]
    to_remove <- (p > alpha) | (abspcor < min_abspcor)
    if (any(to_remove)) {
      rem <- idx[to_remove, , drop = FALSE]
      mirror_rem <- rem[, 2:1, drop = FALSE]
      G[rem] <- G[mirror_rem] <- FALSE
      if (!is.null(sepSet)) {
        pos <- which(needs_update)[to_remove]
        pb <- pbmcapply::progressBar(min = 0, max = length(pos))
        for (k in seq_along(pos)) {
          setTxtProgressBar(pb, value = k)
          a <- rem[k, 1]
          b <- rem[k, 2]
          sep <- res[[pos[k]]]$sepset
          sepSet[[a]][[b]] <- sepSet[[b]][[a]] <- sep
        }
        close(pb)
      }
    }
  }
  done <- !any(needs_update)
  return(list(
    G = G, pMax = pMax, chisqMin = chisqMin,
    absPcorMin = absPcorMin, Threshold = Threshold,
    sampleSize = sampleSize, sepSet = sepSet, done = done
  ))
}


#' Infer graph skeleton via conditional independence tests
#'
#' Constructs graph skeleton by performing iterative conditional independence
#' (CI) tests on gene expression data.
#' Edges are removed based on statistical independence.
#'
#' @param count scRNA-seq count matrix (cells x genes).
#' @param Y scRNA-seq normalized expression matrix (cells × genes). For
#' example, log1p(total UMI corrected count).
#' @param alpha Significance level for CI tests.
#' @param min_abspcor Minimum absolute value of partial correlation for kept edges.
#' @param ncores Number of CPU cores for parallel processing.
#' @param G Optional initial adjacency matrix (defaults to a fully connected
#' graph without self-loops).
#' @param max_order Maximum conditioning set size (0 or 1, default is 1).
#' @param max_thr Maximum threshold for conditional variable (default is 10).
#' @param min_n1 Minimum number of samples satisfying Yk > selected threshold (default is 1000).
#' @param min_n2 Minimum number of samples satisfying Yk > selected threshold,
#' Yi > 0, and Yj > 0 (default is 200).
#' @param sepset Return separation set or not (default is \code{TRUE}).
#'
#' @return A list with:
#' \describe{
#'   \item{graph}{An igraph object representing the inferred skeleton with edge
#'   attributes `pMax`, `chisqMin`, and `absPcorMin`.}
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
#' rownames(count) <- paste0("cell", seq_len(nrow(count)))
#' Y <- log1p(count)
#' skel <- infer_skeleton(
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
infer_skeleton <- function(
    count, Y, alpha, min_abspcor, ncores, G = NULL,
    max_order = 1, max_thr = 10, min_n1 = 1000, min_n2 = 200,
    sepset = TRUE
) {
  .check_skeleton_data(count = count, Y = Y, G = G)
  .check_skeleton_params(
    alpha = alpha,
    min_abspcor = min_abspcor,
    max_order = max_order,
    max_thr = max_thr,
    min_n1 = min_n1,
    min_n2 = min_n2,
    sepset = sepset
  )
  ncores <- .check_ncores(ncores)

  genes <- colnames(Y)
  if (is.null(G)) {
    G <- matrix(TRUE, length(genes), length(genes), dimnames = list(genes, genes))
    diag(G) <- FALSE
  } else {
    G <- G != 0
  }
  last_nedge <- nedge <- sum(G) / 2
  message(paste0("Number of edges in initial graph: ", nedge))
  # Initialize CI test statistics
  pMax <- chisqMin <- absPcorMin <- Threshold <- sampleSize <- matrix(
    NA, length(genes), length(genes), dimnames = list(genes, genes)
  )
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
    res <- .perform_ci_test(
      G = G,
      order = order,
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
      Threshold = Threshold,
      sampleSize = sampleSize,
      sepSet = sepSet,
      ncores = ncores
    )
    G <- res$G
    pMax <- res$pMax
    chisqMin <- res$chisqMin
    absPcorMin <- res$absPcorMin
    Threshold <- res$Threshold
    sampleSize <- res$sampleSize
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
    G = G, pMax = pMax, chisqMin = chisqMin,
    absPcorMin = absPcorMin, Threshold = Threshold,
    sampleSize = sampleSize
  )
  return(list(graph = graph, sepSet = sepSet))
}
