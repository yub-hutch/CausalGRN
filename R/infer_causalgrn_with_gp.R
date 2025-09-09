#' Direct edges using perturbation effect with gene program
#'
#' Infers causal GRN by directing edges based on differential expression caused by perturbations.
#'
#' @param graph Initial igraph object.
#' @param stat Perturbation effect from \code{\link{calc_perturbation_effect}}.
#' @param alpha Numeric representing DE Q-value threshold.
#' @param pname Gene program name.
#' @param pgenes Gene program genes.
#' @param allow_pgene_out Allow program genes as parents of non-program genes or not.
#' @param conservative Logical indicating whether to make conservative inference (Default is \code{TRUE}).
#'
#' @return igraph object.
#' @export
infer_causalgrn_with_gp <- function(graph, stat, alpha, pname, pgenes, allow_pgene_out, conservative = TRUE) {
  kos <- unique(stat$ko)
  genes <- unique(stat$gene)
  stopifnot(all(kos %in% genes))
  stopifnot(nrow(stat) == length(kos) * length(genes))
  stopifnot(setequal(genes, igraph::V(graph)$name))
  # Extract DE adjusted p-values
  adj_pv_mat <- as.matrix(stats::xtabs(adj_pv ~ ko + gene, data = stat))
  # Order 1 orientation
  edges_to_delete <- c()
  visited_kos <- c()
  for (ko in kos) {
    nodes <- igraph::neighbors(graph, ko, mode = 'all')$name
    nodes <- setdiff(nodes, visited_kos)
    if (ko %in% pgenes) nodes <- setdiff(nodes, pname)
    visited_kos <- c(visited_kos, ko)
    if (length(nodes) == 0) next
    for (node in nodes) {
      adj_pv <- adj_pv_mat[ko, node]
      if (node %in% kos) {
        adj_pv_reverse <- adj_pv_mat[node, ko]
        if (adj_pv < alpha && adj_pv_reverse > alpha) {
          edges_to_delete <- c(edges_to_delete, c(node, ko))
        } else if (adj_pv > alpha && adj_pv_reverse < alpha) {
          edges_to_delete <- c(edges_to_delete, c(ko, node))
        }
      } else {
        if (adj_pv < alpha) {
          edges_to_delete <- c(edges_to_delete, c(node, ko))
        } else if (!conservative) {
          edges_to_delete <- c(edges_to_delete, c(ko, node))
        }
      }
    }
  }
  # If not allow program genes as parents of non program genes
  if (!allow_pgene_out) {
    npgenes <- setdiff(genes, c(pname, pgenes))
    if (length(npgenes)) {
      from <- rep(pgenes, each = length(npgenes))
      to <- rep(npgenes, times = length(pgenes))
      edges_to_delete = c(edges_to_delete, c(rbind(from, to)))
    }
  }
  eids_to_delete <- setdiff(igraph::get_edge_ids(graph, edges_to_delete, directed = TRUE), 0)
  if (length(eids_to_delete)) graph <- igraph::delete_edges(graph, eids_to_delete)
  return(graph)
}
