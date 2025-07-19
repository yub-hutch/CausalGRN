#' Calculate Perturbation Effect
#'
#' Calculate differential expression statistics between wild-type cells and perturbed cells for each perturbation.
#'
#' @param Y Matrix of normalized scRNA-seq data of wild-type and perturbed cells.
#' @param group Named vector of cell label: 'WT' for wild-type cells, perturbed gene for perturbed cell.
#' @param ncores Number of CPUs to use.
#' @return Tibble with columns:
#' \itemize{
#'   \item \code{ko}: Perturbed gene.
#'   \item \code{gene}: Affected gene.
#'   \item \code{diff}: Mean difference (Perturbed - WT).
#'   \item \code{cd}: Cohen's D (Normalized mean difference).
#'   \item \code{wilcox_pv}: Wilcoxon rank sum test P-value.
#'   \item \code{t_pv}: T-test P-value.
#'   \item \code{cor_pearson}: Pearson correlation in wild-type cells.
#'   \item \code{cor_spearman}: Spearman correlation in all cells.
#'   \item \code{adj_wilcox_pv}: BH-adjusted Wilcoxon rank sum test P-value.
#'   \item \code{adj_t_pv}: BH-adjusted T-test P-value.
#'   \item \code{adj_pv}: The larger of `adj_wilcox_pv` and `adj_wilcox_pv`.
#' }
#'
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
#' print(stat)
#' @export
calc_perturbation_effect <- function(Y, group, ncores) {
  # Check group
  stopifnot(identical(rownames(Y), names(group)))
  stopifnot('group must contain wild-type cells labeled as WT' = 'WT' %in% group)
  stopifnot('Remove groups with less than 50 cells' = all(table(group) >= 50))
  # Check potential memory issue
  kos <- setdiff(group, 'WT')
  genes <- colnames(Y)
  use_batches <- nrow(Y) > 10e3 && ncol(Y) > 1e3 && ncores >= 10
  if (!use_batches) {
    wt <- Y[group == 'WT', ]
    stat <- do.call(rbind, pbmcapply::pbmcmapply(\(ko) {
      pt <- Y[group == ko, ]
      diffs <- colMeans(pt) - colMeans(wt)
      pooled_sds <- apply(rbind(wt, pt), 2, sd)
      cds <- diffs / pooled_sds
      wilcox_pvs <- matrixTests::col_wilcoxon_twosample(wt, pt, exact = F, correct = T)$pvalue
      t_pvs <- matrixTests::col_t_welch(wt, pt)$pvalue
      cors_pearson <- c(cor(wt[, ko], wt))
      cors_spearman <- c(cor(wt[, ko], wt, method = 'spearman'))
      dplyr::tibble(
        ko = ko, gene = genes,
        diff = diffs, cd = cds,
        wilcox_pv = wilcox_pvs, t_pv = t_pvs,
        cor_pearson = cors_pearson, cor_spearman = cors_spearman
      )
    }, kos, SIMPLIFY = F, mc.cores = ncores))
  } else {
    message('Calculating DE statistics in batches ...')
    ko_expr_wt_list <- sapply(kos, simplify = FALSE, \(ko) Y[group == 'WT', ko])
    chunks <- split(seq_along(genes), f = ceiling(seq_along(genes) / 200))
    message(paste0('Number of batches: ', length(chunks)))
    message('------------------------------------------------------')
    stat <- do.call(rbind, lapply(seq_along(chunks), \(i) {
      message(paste0('Batch ', i))
      cols <- chunks[[i]]
      sub_Y <- Y[, cols, drop = FALSE] # Read data from disk to memory
      wt <- sub_Y[group == 'WT', , drop = FALSE]
      do.call(rbind, pbmcapply::pbmcmapply(\(ko, ko_expr_wt) {
        pt <- sub_Y[group == ko, , drop = FALSE]
        diffs <- colMeans(pt) - colMeans(wt)
        pooled_sds <- apply(rbind(wt, pt), 2, sd)
        cds <- diffs / pooled_sds
        wilcox_pvs <- matrixTests::col_wilcoxon_twosample(wt, pt, exact = FALSE, correct = TRUE)$pvalue
        t_pvs <- matrixTests::col_t_welch(wt, pt)$pvalue
        cors_pearson <- c(cor(ko_expr_wt, wt))
        cors_spearman <- c(cor(ko_expr_wt, wt, method = 'spearman'))
        dplyr::tibble(
          ko = ko, gene = genes[cols],
          diff = diffs, cd = cds,
          wilcox_pv = wilcox_pvs, t_pv = t_pvs,
          cor_pearson = cors_pearson, cor_spearman = cors_spearman
        )
      }, kos, ko_expr_wt_list, SIMPLIFY = FALSE, mc.cores = ncores))
    }))
  }
  stat$wilcox_adj_pv = p.adjust(stat$wilcox_pv, method = 'BH')
  stat$t_adj_pv = p.adjust(stat$t_pv, method = 'BH')
  stat$adj_pv = pmax(stat$wilcox_adj_pv, stat$t_adj_pv)
  return(stat)
}
