.calc_line <- function(cors, cds, cor_probs) {
  median_cd_probs <- sapply(cor_probs, function(cor_prob) {
    cor_thr <- quantile(abs(cors), probs = cor_prob)
    curr_cds <- abs(cds)[abs(cors) >= cor_thr]
    median_cd <- median(curr_cds)
    return(mean(median_cd >= abs(cds)))
  })
  return(tibble::tibble(cor_prob = cor_probs, median_cd_prob = median_cd_probs))
}

.calc_auc <- function(line) {
  line <- line[order(line$cor_prob), ]
  x <- line$cor_prob
  y <- line$median_cd_prob
  dx <- diff(x)
  midy <- (y[-1] + y[-length(y)]) / 2
  auc <- sum(dx * midy)
  return(auc)
}

#' Calculate AUC measuring consistency between correlation and differential expression
#'
#' Calculate AUC measuring consistency between Pearson/Spearman correlations and Cohen' D, under specified threshold.
#'
#' @param stat DE statistics, see \code{\link{calc_perturbation_effect}}.
#' @param thr Numeric value in [0, 1) indicating the minimum quantile of absolute correlation to include in AUC calculation.
#' @param ncores Number of cores (default is 1).
#'
#' @return AUC matrix with rows representing KOs and two AUC columns based on 'pearson' and 'spearman' correlations.
#' @examples
#' # x -> y -> z
#' set.seed(123)
#' n <- 1e3
#' x0 <- rnorm(n)
#' x1 <- rnorm(n, mean = -1)
#' x <- c(x0, x1)
#' group <- c(rep('WT', n), rep('x', n))
#' y <- rnorm(2 * n, mean = x, sd = 1)
#' z <- rnorm(2 * n, mean = y, sd = 1)
#' data <- cbind(x = x, y = y, z = z)
#' stat <- calc_perturbation_effect(
#'   Y = data,
#'   group = group,
#'   ncores = 2
#' )
#' calc_auc(stat, thr = 0)
#' @export
calc_auc <- function(stat, thr, ncores = 1) {
  stopifnot(thr >= 0 & thr < 1)
  kos <- unique(stat$ko)
  cor_probs <- seq(thr, 1, by = 0.01)
  auc <- do.call(rbind, pbmcapply::pbmclapply(kos, function(ko) {
    cds <- stat$cd[(stat$ko == ko) & (stat$gene != ko)]
    cors_pearson <- stat$cor_pearson[(stat$ko == ko) & (stat$gene != ko)]
    cors_spearman <- stat$cor_spearman[(stat$ko == ko) & (stat$gene != ko)]
    line_pearson <- .calc_line(cors_pearson, cds, cor_probs)
    line_spearman <- .calc_line(cors_spearman, cds, cor_probs)
    auc_pearson <- .calc_auc(line_pearson) / (1 - thr)
    auc_spearman <- .calc_auc(line_spearman) / (1 - thr)
    c(pearson = auc_pearson, spearman = auc_spearman)
  }, mc.cores = ncores))
  rownames(auc) <- kos
  return(auc)
}

