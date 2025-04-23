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
#'   \item \code{cor_wt}: Correlation in wild-type cells.
#'   \item \code{cor_all}: Correlation in all cells.
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
  use_disk <- nrow(Y) > 10e3 && ncol(Y) > 1e3 && ncores >= 10
  kos <- setdiff(group, 'WT')
  genes <- colnames(Y)
  if (!use_disk) {
    wt <- Y[group == 'WT', ]
    stat <- do.call(rbind, pbmcapply::pbmcmapply(function(ko) {
      pt <- Y[group == ko, ]
      diffs <- colMeans(pt) - colMeans(wt)
      pooled_sds <- apply(rbind(wt, pt), 2, sd)
      cds <- diffs / pooled_sds
      wilcox_pvs <- matrixTests::col_wilcoxon_twosample(wt, pt, exact = F, correct = T)$pvalue
      t_pvs <- matrixTests::col_t_welch(wt, pt)$pvalue
      cors_wt <- c(cor(wt[, ko], wt))
      cors_all <- c(cor(Y[, ko], Y))
      dplyr::tibble(
        ko = ko, gene = genes,
        diff = diffs, cd = cds,
        wilcox_pv = wilcox_pvs, t_pv = t_pvs,
        cor_wt = cors_wt, cor_all = cors_all
      )
    }, kos, SIMPLIFY = F, mc.cores = ncores))
  } else {
    message('Y is a big matrix. Writing to disk ...')
    ngene <- ncol(Y)
    genes <- colnames(Y)
    Y <- bigstatsr::as_FBM(Y, backingfile = 'Y_fbm') # Row and column names erased
    gc()
    message('Calculating DE statistics ...')
    # Do by gene chunks
    chunk_size <- 200
    chunks <- split(seq_len(ngene), f = ceiling(seq_len(ngene) / chunk_size))
    pb <- pbmcapply::progressBar(min = 0, max = length(chunks))
    stat <- do.call(rbind, lapply(seq_along(chunks), function(i) {
      setTxtProgressBar(pb = pb, value = i)
      cols <- chunks[[i]]
      sub_Y <- Y[, cols, drop = FALSE] # Read data from disk to memory
      wt <- sub_Y[group == 'WT', , drop = FALSE]
      sub_stat <- do.call(rbind, parallel::mclapply(kos, function(ko) {
        ko_col <- match(ko, genes)
        ko_expr_wt <- Y[group == 'WT', ko_col]
        ko_expr_all <- Y[, ko_col]
        pt <- sub_Y[group == ko, , drop = FALSE]
        diffs <- colMeans(pt) - colMeans(wt)
        pooled_sds <- apply(rbind(wt, pt), 2, sd)
        cds <- diffs / pooled_sds
        wilcox_pvs <- matrixTests::col_wilcoxon_twosample(wt, pt, exact = F, correct = T)$pvalue
        t_pvs <- matrixTests::col_t_welch(wt, pt)$pvalue
        cors_wt <- c(cor(ko_expr_wt, wt))
        cors_all <- c(cor(ko_expr_all, sub_Y))
        dplyr::tibble(
          ko = ko, gene = genes[cols],
          diff = diffs, cd = cds,
          wilcox_pv = wilcox_pvs, t_pv = t_pvs,
          cor_wt = cors_wt, cor_all = cors_all
        )
      }, mc.cores = ncores))
      return(sub_stat)
    }))
    close(pb)
  }
  stat$wilcox_adj_pv = p.adjust(stat$wilcox_pv, method = 'BH')
  stat$t_adj_pv = p.adjust(stat$t_pv, method = 'BH')
  stat$adj_pv = pmax(stat$wilcox_adj_pv, stat$t_adj_pv)
  return(stat)
}
