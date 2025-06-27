#' Evaluate Structure Difference to True Graph
#'
#' This function evaluates the structure difference between the estimated directed graph and the true graph.
#'
#' @param g igraph object of the estimated graph.
#' @param g0 igraph object of the true graph.
#' @return Tibble with columns:
#' \itemize{
#'   \item \code{nedge_true}: Number of edges in the true graph.
#'   \item \code{nedge}: Number of edges in the estimated graph.
#'   \item \code{intersection}: Number of overlap edges.
#'   \item \code{jaccard}: Jaccard index between the edges of the two graphs.
#'   \item \code{recall}: Recall of the estimated graph.
#'   \item \code{precision}: Precision of the estimated graph.
#'   \item \code{f1}: F1 score of the estimated graph.
#'   \item \code{shd}: Structural Hamming Distance (SHD) between the true and estimated graphs.
#' }
#' @export
structure_difference <- function(g, g0) {
  dplyr::tibble(
    nedge_true = igraph::ecount(g0),
    nedge = igraph::ecount(g),
    intersection = igraph::ecount(igraph::intersection(g0, g)),
    jaccard = intersection / igraph::ecount(igraph::union(g0, g)),
    recall = intersection / nedge_true,
    precision = intersection / nedge,
    f1 = 2 * precision * recall / (precision + recall),
    shd = nedge_true + nedge - 2 * intersection
  )
}


#' Evaluate Structure Difference to True Graph at Different Breaks of Top Predicted Edges
#'
#' This function extract the subgraph with top predicted edges and compare its structure with the true graph.
#'
#' @param g0 igraph object of the true graph.
#' @param pred Data frame of predicted directed edges with column \code{score} representing edge ranking.
#' @param ntops Vector of number of top predicted edges to evaluate.
#' @return Tibble with columns:
#' \itemize{
#'   \item \code{ntop}: Number of top predicted edges evaluated.
#'   \item \code{nedge_true}: Number of edges in the true graph.
#'   \item \code{nedge}: Number of edges in the estimated graph.
#'   \item \code{intersection}: Number of overlap edges.
#'   \item \code{jaccard}: Jaccard index between the edges of the two graphs.
#'   \item \code{recall}: Recall of the estimated graph.
#'   \item \code{precision}: Precision of the estimated graph.
#'   \item \code{f1}: F1 score of the estimated graph.
#'   \item \code{shd}: Structural Hamming Distance (SHD) between the true and estimated graphs.
#' }
#' @export
structure_difference_at_breaks <- function(pred, g0, ntops) {
  stopifnot('score' %in% names(pred))
  pred = pred %>% dplyr::arrange(dplyr::desc(score))
  ntops = sort(ntops[ntops <= nrow(pred)])
  do.call(rbind, lapply(ntops, function(ntop) {
    g = igraph::graph_from_data_frame(pred[seq(ntop), ])
    metrics = structure_difference(g = g, g0 = g0)
    dplyr::tibble(ntop = ntop, metrics)
  }))
}
