#' Scatter plot of Correlation versus Cohens' D
#'
#' Draw scatter plot between Cohen's D and correlation in wild-type cells of KO-gene pairs.
#'
#' @param stat Calculated from \code{\link{calc_perturbation_effect}}.
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
plot_cohens_d_vs_cor <- function(stat) {
  plot_data <- stat |> dplyr::filter(ko != gene)
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(cor_wt, abs(cd))) +
    ggpointdensity::geom_pointdensity() +
    ggplot2::scale_color_viridis_c(guide = "none") +
    ggplot2::labs(x = 'Correlation in WT cells', y = "|Cohens' D|") +
    ggplot2::geom_hline(yintercept = median(abs(plot_data$cd)), color = 'red', linetype = 2) +
    ggplot2::geom_vline(xintercept = 0, color = 'red', linetype = 2) +
    ggpubr::theme_pubr()
  return(p)
}
