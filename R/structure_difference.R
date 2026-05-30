#' Evaluate Structure Difference to True Graph
#'
#' This function evaluates the structure difference between the estimated
#' directed graph and the true graph.
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
#'   \item \code{shd}: Structural Hamming Distance (SHD), as computed by
#'     \code{pcalg::shd()}.
#' }
#' @export
structure_difference <- function(g, g0) {
  if (!inherits(g, "igraph") || !inherits(g0, "igraph")) {
    stop("'g' and 'g0' must be igraph objects.", call. = FALSE)
  }

  nodes <- igraph::V(g0)$name
  if (is.null(nodes) || anyNA(nodes) || any(nodes == "")) {
    stop("'g0' must have non-missing vertex names.", call. = FALSE)
  }

  true_adj <- igraph::as_adjacency_matrix(g0, sparse = FALSE) > 0
  true_adj <- true_adj[nodes, nodes, drop = FALSE]

  pred_adj <- matrix(
    FALSE,
    length(nodes),
    length(nodes),
    dimnames = list(nodes, nodes)
  )
  pred_edges <- matrix(character(), ncol = 2)
  if (igraph::ecount(g) > 0L) {
    pred_edges <- igraph::as_edgelist(g, names = TRUE)
    if (!all(pred_edges %in% nodes)) {
      stop("'g' edges must only use vertices in 'g0'.", call. = FALSE)
    }
    pred_adj[pred_edges] <- TRUE
  }

  nedge_true <- sum(true_adj)
  nedge <- sum(pred_adj)
  intersection <- sum(true_adj & pred_adj)
  union_edges <- sum(true_adj | pred_adj)
  jaccard <- intersection / union_edges
  recall <- intersection / nedge_true
  precision <- intersection / nedge
  f1 <- 2 * intersection / (nedge_true + nedge)

  pred_edge_df <- unique(data.frame(
    from = pred_edges[, 1],
    to = pred_edges[, 2]
  ))
  pred_graph <- igraph::graph_from_data_frame(
    pred_edge_df,
    directed = TRUE,
    vertices = data.frame(name = nodes)
  )
  shd <- pcalg::shd(
    igraph::as_graphnel(pred_graph),
    igraph::as_graphnel(g0)
  )

  dplyr::tibble(
    nedge_true = nedge_true,
    nedge = nedge,
    intersection = intersection,
    jaccard = jaccard,
    recall = recall,
    precision = precision,
    f1 = f1,
    shd = shd
  )
}


#' Evaluate Structure Difference at Top-Edge Breaks
#'
#' This function extracts the subgraph with top predicted edges and compares its
#' structure with the true graph.
#'
#' @param g0 igraph object of the true graph.
#' @param pred Data frame of predicted directed edges with columns \code{from},
#'   \code{to}, and \code{score}.
#' @param ntops Vector of number of top predicted edges to evaluate.
#' @param ncores Number of cores to use for parallel computation (default is 1).
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
#'   \item \code{shd}: Structural Hamming Distance (SHD), as computed by
#'     \code{pcalg::shd()}.
#' }
#' @export
structure_difference_at_breaks <- function(pred, g0, ntops, ncores = 1) {
  required_cols <- c("from", "to", "score")
  if (!is.data.frame(pred) || !all(required_cols %in% names(pred))) {
    stop(
      "'pred' must contain columns 'from', 'to', and 'score'.",
      call. = FALSE
    )
  }
  .check_ncores(ncores)

  pred <- pred[
    order(pred$score, decreasing = TRUE),
    required_cols,
    drop = FALSE
  ]
  ntops <- sort(ntops[ntops <= nrow(pred)])

  metrics_list <- .parallel_lapply(
    ntops,
    function(ntop) {
      g <- igraph::graph_from_data_frame(
        pred[seq_len(ntop), , drop = FALSE],
        directed = TRUE,
        vertices = data.frame(name = igraph::V(g0)$name)
      )
      metrics <- structure_difference(g = g, g0 = g0)
      dplyr::tibble(ntop = ntop, metrics)
    },
    ncores = ncores,
    export = c("pred", "g0", "structure_difference")
  )

  do.call(rbind, metrics_list)
}
