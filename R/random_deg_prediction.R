#' Random‑null correlation for a leave‑one‑group‑out baseline
#'
#' Return the analytic (population‑level) Pearson correlation between the true differential expression
#' (test‑KO – WT) and the trivial baseline that uses the weighted mean of all training cells (WT + remaining KOs) as a predictor.
#'
#' @param group_sizes Integer vector of group (cell) counts. The first element must be WT; the rest are KO conditions.
#' @param test_idx Integer index of the KO group that is held out.
#' @return The analytic random‑null correlation.
#' @examples
#' # Three groups: WT = 1000 cells, KO‑1 = 100, KO‑2 = 100 (held out)
#' group_sizes <- c(1000, 100, 100)
#' test_idx <- 3
#' calc_random_cor(group_sizes, test_idx)
#' @export
calc_random_cor <- function(group_sizes, test_idx) {
  wt_idx <- 1
  n0 <- group_sizes[wt_idx]
  nt <- group_sizes[test_idx]
  ntrain <- group_sizes[-test_idx]
  wtrain <- ntrain / sum(ntrain)
  w0 <- n0 / sum(ntrain)
  A <- 1 / n0 + 1 / nt
  B <- (1 - w0)^2 / n0 + sum(wtrain[-wt_idx]^2 / group_sizes[-c(wt_idx, test_idx)])
  return((1 - w0) / n0 / sqrt(A * B))
}


#' Monte‑Carlo check of the random‑null correlation
#'
#' Draw full cell × gene matrices from the global‑mean‑centered null model, repeats \code{n_runs} independent experiments,
#' and returns the vector of sample correlations. The empirical mean should hover around \code{\link{calc_random_cor}}.
#'
#' @param group_sizes Integer vector of cell counts (WT first).
#' @param test_idx    Integer index of the held‑out KO.
#' @param n_genes     Number of genes to simulate in every run.
#' @param n_runs      Number of independent Monte‑Carlo experiments.
#' @param ncores      Number of CPU cores.
#' @return Numeric vector of length `n_runs` with sample correlations.
#' @examples
#' group_sizes <- c(1000, 100, 100)
#' test_idx <- 3
#' r0 <- calc_random_cor(group_sizes, test_idx)
#' emp_r0s <- simulate_random_cor(group_sizes, test_idx, n_genes = 1000, n_runs = 100, ncores = 4)
#' cat(sprintf("ρ₀ analytic = %.3f\n", r0))
#' cat(sprintf("Empirical r = %.3f  (sd = %.3f)\n", mean(emp_r0s), sd(emp_r0s)))
#' @export
simulate_random_cor <- function(group_sizes, test_idx, n_genes, n_runs, ncores) {
  ngroup <- length(group_sizes)
  group <- rep(seq_len(ngroup), times = group_sizes)
  cors <- unlist(pbmcapply::pbmclapply(seq_len(n_runs), function(k) {
    expr <- matrix(rnorm(sum(group_sizes) * n_genes), nrow = sum(group_sizes), ncol = n_genes)
    expr <- sweep(expr, 2, colMeans(expr))
    mu_wt <- colMeans(expr[group == 1, , drop = FALSE])
    mu_t <- colMeans(expr[group == test_idx, , drop = FALSE])
    train_cells <- group != test_idx
    mu_train <- colMeans(expr[train_cells, , drop = FALSE])
    de <- mu_t - mu_wt
    pred_de <- mu_train - mu_wt
    return(cor(de, pred_de))
  }, mc.cores = ncores))
  return(cors)
}
