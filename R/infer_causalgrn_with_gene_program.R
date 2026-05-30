#' Direct edges using perturbation effect with a gene program node
#'
#' Infers causal GRN by directing edges based on differential expression caused
#' by perturbations. Program genes are always allowed as parents of
#' non-program genes.
#'
#' @param graph Initial igraph object.
#' @param stat Perturbation effect from \code{\link{calc_perturbation_effect}}.
#' @param alpha Numeric representing DE Q-value threshold.
#' @param pname Gene program node name.
#' @param pgenes Gene program member genes.
#' @param conservative Logical indicating whether to make conservative inference
#' (Default is \code{TRUE}).
#' @param max_order Integer representing the maximum order for DE descendant
#' inference. Can be 1 or 2. Default is 1.
#'
#' @return igraph object.
#' @export
infer_causalgrn_with_gene_program <- function(
    graph, stat, alpha, pname, pgenes, conservative = TRUE, max_order = 1
) {
  .check_causalgrn_data(graph = graph, stat = stat)
  .check_causalgrn_params(
    alpha = alpha,
    conservative = conservative,
    max_order = max_order
  )

  nodes <- igraph::V(graph)$name
  genes <- setdiff(nodes, pname)
  stopifnot(
    length(pname) == 1L,
    pname %in% nodes,
    is.character(pgenes),
    all(pgenes %in% genes)
  )

  kos <- unique(stat$ko)

  # Extract DE adjusted p-values
  adj_pv_mat <- as.matrix(stats::xtabs(adj_pv ~ ko + gene, data = stat))

  # Order 1 orientation
  edges_to_delete <- character(0)
  visited_kos <- character(0)
  for (ko in kos) {
    nodes <- igraph::neighbors(graph, ko, mode = 'all')$name
    nodes <- setdiff(nodes, visited_kos)
    if (ko %in% pgenes) {
      nodes <- setdiff(nodes, pname)
    }
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
  edge_ids_to_delete <- setdiff(
    igraph::get_edge_ids(graph, edges_to_delete, directed = TRUE),
    0
  )
  if (length(edge_ids_to_delete)) {
    graph <- igraph::delete_edges(graph, edge_ids_to_delete)
  }

  # Order 2 orientation
  if (max_order == 2) {
    edges_to_delete <- character(0)
    for (ko in kos) {
      children <- setdiff(
        igraph::neighbors(graph, ko, mode = 'out')$name,
        igraph::neighbors(graph, ko, mode = 'in')$name
      )
      children <- setdiff(children, kos)
      if (length(children) == 0) next
      for (child in children) {
        order2_nodes <- intersect(
          igraph::neighbors(graph, child, mode = 'in')$name,
          igraph::neighbors(graph, child, mode = 'out')$name
        )
        if (length(order2_nodes) == 0) next
        for (order2_node in order2_nodes) {
          order2_adj_pv <- adj_pv_mat[ko, order2_node]
          if (order2_adj_pv < alpha) {
            distance_wo_child <- igraph::distances(
              graph = igraph::delete_vertices(graph, child),
              v = ko,
              to = order2_node,
              mode = 'out'
            )[1, 1]
            is_child_on_every_path <- is.infinite(distance_wo_child)
            if (is_child_on_every_path) {
              edges_to_delete <- c(edges_to_delete, paste0(order2_node, '->', child))
            }
          }
        }
      }
    }

    # Deduplicate edge deletion candidates.
    if (length(edges_to_delete)) {
      edges_to_delete <- unique(edges_to_delete)
    }

    # Drop edges with conflict
    if (length(edges_to_delete)) {
      rev_edges_to_delete <- sub("^(.*)->(.*)$", "\\2->\\1", edges_to_delete)
      to_drop <- rev_edges_to_delete %in% edges_to_delete
      edges_to_delete <- edges_to_delete[!to_drop]
    }
    # Delete remaining edges
    if (length(edges_to_delete)) {
      edges_to_delete <- unlist(strsplit(edges_to_delete, '->'))
      edge_ids_to_delete <- igraph::get_edge_ids(
        graph,
        edges_to_delete,
        directed = TRUE
      )
      graph <- igraph::delete_edges(graph, edge_ids_to_delete)
    }
  }

  return(graph)
}
