#' Scatter plot of Correlation versus Cohens' D
#'
#' Draw scatter plot between Cohen's D and correlation in wild-type cells of KO-gene pairs.
#'
#' @param stat Calculated from \code{\link{calc_perturbation_effect}}.
#' @param cor Correlation across "WT" or "All" cells.
#' @return ggplot object.
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
#' plot_cohens_d_vs_cor(stat)
#'
#' @export
plot_cohens_d_vs_cor <- function(stat, cor = c('WT', 'All')) {
  cor <- match.arg(cor)
  if (cor == 'WT') {
    plot_data <- stat |> dplyr::filter(ko != gene) |> dplyr::mutate(cor = cor_wt)
    xlab <- 'Correlation in WT cells'
  } else {
    plot_data <- stat |> dplyr::filter(ko != gene) |> dplyr::mutate(cor = cor_all)
    xlab <- 'Correlation in all cells'
  }
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(cor, abs(cd))) +
    ggpointdensity::geom_pointdensity() +
    ggplot2::scale_color_viridis_c(guide = "none") +
    ggplot2::labs(x = xlab, y = "|Cohens' D|") +
    ggplot2::geom_hline(yintercept = median(abs(plot_data$cd)), color = 'red', linetype = 2) +
    ggplot2::geom_vline(xintercept = 0, color = 'red', linetype = 2) +
    ggpubr::theme_pubr()
  return(p)
}


#' Plot adjusted P-values for KO gene pairs in -log10 scale
#'
#' @param stat Tibble containing DE statistics (See \code{\link{calc_perturbation_effect}}).
#' @param eps Numeric, minimum value to clamp p-values to. Default is 1e-5.
#' @return ggplot object.
#' @export
plot_adj_pv_for_ko_pairs <- function(stat, eps = 1e-5) {
  stat1 <- stat |> dplyr::filter(ko > gene)
  stat2 <- stat |> dplyr::filter(ko < gene)
  plot_data <- stat1 |>
    dplyr::inner_join(stat2, by = c(ko = "gene", gene = "ko"), suffix = c("_1", "_2")) |>
    dplyr::mutate(
      adj_pv_1_clamped = pmax(adj_pv_1, eps),
      adj_pv_2_clamped = pmax(adj_pv_2, eps)
    )
  thr <- -log10(0.05)
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(-log10(adj_pv_1_clamped), -log10(adj_pv_2_clamped))) +
    ggpointdensity::geom_pointdensity() +
    ggplot2::scale_color_viridis_c(guide = "none") +
    ggplot2::labs(
      x = "-log10 adjusted p-value of gene A when gene B is perturbed",
      y = "-log10 adjusted p-value of gene B when gene A is perturbed"
    ) +
    ggplot2::geom_vline(xintercept = thr, color = "red", linetype = 2) +
    ggplot2::geom_hline(yintercept = thr, color = "red", linetype = 2) +
    ggpubr::theme_pubr()
  return(p)
}
