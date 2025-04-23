#' Direct edges using perturbation effect
#'
#' Infers causal GRN by directing edges based on differential expression caused by perturbations.
#'
#' @param graph Initial igraph object.
#' @param stat Perturbation effect from \code{\link{calc_perturbation_effect}}.
#' @param alpha Numeric representing DE Q-value threshold.
#' @param conservative Logical indicating whether to make conservative inference (Default is \code{TRUE}).
#' @param max_order Integer representing the maximum order for DE descendant inference. Can be 1 or 2. Default is 1.
#' @param max_dist Integer representing the maximum distance allowed for a perturbed gene to cause DE on another gene. Ignore for `max_order = 1`. Default is Inf.
#' @param evidence Integer representing number of evidences needed for second-order orientation. Ignore for `max_order = 1`. Default is 1.
#'
#' @examples
#' # No hub genes
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
#' skel <- infer_skeleton_from_perturbed(
#'   Y = data,
#'   group = group,
#'   alpha = 1e-3,
#'   ncores = 2,
#'   max_order = 1
#' )
#' stat <- calc_perturbation_effect(
#'   Y = data,
#'   group = group,
#'   ncores = 2
#' )
#' cgrn <- infer_causalgrn(
#'   graph = skel$graph,
#'   stat = stat,
#'   alpha = 0.05,
#'   max_order = 2
#' )
#' plot(cgrn)
#' print(igraph::as_data_frame(cgrn))
#'
#' @return igraph object.
#' @export
infer_causalgrn <- function(graph, stat, alpha, conservative = TRUE, max_order = 1, max_dist = Inf, evidence = 1) {
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
    visited_kos <- c(visited_kos, ko)
    if (length(nodes) == 0) {
      next
    }
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
  edge_ids_to_delete <- igraph::get.edge.ids(graph, edges_to_delete, directed = TRUE)
  edge_ids_to_delete <- edge_ids_to_delete[edge_ids_to_delete != 0]
  if (length(edge_ids_to_delete) > 0) {
    graph <- igraph::delete_edges(graph, edge_ids_to_delete)
  }
  # Order 2 orientation
  if (max_order == 2) {
    edges_to_delete <- c()
    for (ko in kos) {
      children <- setdiff(
        igraph::neighbors(graph, ko, mode = 'out')$name,
        igraph::neighbors(graph, ko, mode = 'in')$name
      )
      children <- setdiff(children, kos)
      if (length(children) == 0) {
        next
      }
      for (child in children) {
        order2_nodes <- intersect(
          igraph::neighbors(graph, child, mode = 'in')$name,
          igraph::neighbors(graph, child, mode = 'out')$name
        )
        if (length(order2_nodes) == 0) {
          next
        }
        for (order2_node in order2_nodes) {
          order2_adj_pv <- adj_pv_mat[ko, order2_node]
          if (order2_adj_pv < alpha) {
            distance_wo_child <- igraph::distances(
              graph = igraph::delete_vertices(graph, child),
              v = ko,
              to = order2_node,
              mode = 'out'
            )[1, 1]
            is_child_on_all_paths <- (distance_wo_child >= max_dist + 1)
            if (is_child_on_all_paths) {
              edges_to_delete <- c(edges_to_delete, paste0(order2_node, '->', child))
            }
          }
        }
      }
    }
    # Only exclude edges with required evidence
    if (length(edges_to_delete) > 0) {
      edges_to_delete <- names(which(table(edges_to_delete) >= evidence))
    }
    # Drop edges with conflict
    if (length(edges_to_delete)) {
      rev_edges_to_delete <- sub("^(.*)->(.*)$", "\\2->\\1", edges_to_delete)
      to_drop <- rev_edges_to_delete %in% edges_to_delete
      edges_to_delete <- edges_to_delete[!to_drop]
    }
    # Delete remaining edges
    if (length(edges_to_delete) > 0) {
      edges_to_delete <- unlist(strsplit(edges_to_delete, '->'))
      edge_ids_to_delete <- igraph::get.edge.ids(graph, edges_to_delete, directed = TRUE)
      graph <- igraph::delete_edges(graph, edge_ids_to_delete)
    }
  }
  return(graph)
}
