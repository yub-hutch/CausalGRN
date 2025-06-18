#' Calculate partial correlation for sparse count data
#'
#' Correct partial correlation for sparse count data by thresholding conditional variable.
#'
#' @param i Integer.
#' @param j Integer.
#' @param k Integer representing index of conditional variable.
#' @param count Count matrix (cells x genes).
#' @param Y Normalized expression matrix (cells x genes).
#' @param max_thr Maximum threshold for conditional variable.
#' @param min_n1 Minimum number of samples satisfying Yk > selected threshold.
#' @param min_n2 Minimum number of samples satisfying Yk > selected threshold, Yi > 0, and Yj > 0.
#' @return Data frame of threshold, sample size, partial correlation estimate, z-score, and p-value.
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
#' # True partial correlation in latent space
#' ppcor::pcor(u)
#' # Uncorrected partial correlation
#' ppcor::pcor.test(Y[, 1], Y[, 3], Y[, 2])
#' ppcor::pcor.test(Y[, 1], Y[, 2], Y[, 3])
#' # Corrected partial correlation
#' calc_pcor(i = 1, j = 3, k = 2, count = count, Y = Y, max_thr = 20, min_n1 = 2000, min_n2 = 400)
#' calc_pcor(i = 1, j = 2, k = 3, count = count, Y = Y, max_thr = 20, min_n1 = 2000, min_n2 = 400)
#' @export
calc_pcor <- function(i, j, k, count, Y, max_thr, min_n1, min_n2) {
  selected_thr <- -1
  for (thr in seq(max_thr, 0)) {
    idx <- count[, k] > thr
    if (sum(idx) < min_n1) next
    if (sum((count[idx, i] > 0) & (count[idx, j] > 0)) < min_n2) next
    selected_thr <- thr
    break
  }
  idx <- count[, k] > selected_thr
  res <- ppcor::pcor.test(x = Y[idx, i], y = Y[idx, j], z = Y[idx, k], method = 'pearson')
  res$threshold <- selected_thr
  return(res)
}
