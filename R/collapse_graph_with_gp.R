#' Collapse graph with gene programs
#'
#' Convert g -> program to g -> all program genes, and drop program.
#'
#' @param graph igraph object.
#' @param pname Gene program name.
#' @param pgenes Gene program genes.
#' @param Y scRNA-seq normalized expression matrix.
#' @param allow_pgene_out Allow program genes as parents of non-program genes or not.
#'
#' @return igraph object.
#' @export
collapse_graph_with_gp <- function(graph, pname, pgenes, Y, allow_pgene_out) {
  genes <- igraph::V(graph)$name
  # Convert g -> program to g -> all program genes
  gs <- igraph::neighbors(graph, pname, mode = 'in')$name
  if (length(gs)) {
    from <- rep(gs, each = length(pgenes))
    to <- rep(pgenes, times = length(gs))
    # Exclude self loop of program genes
    sloop <- from == to
    from <- from[!sloop]
    to <- to[!sloop]
    # For existing edges, modify absPcorMin
    eids <- igraph::get_edge_ids(graph, vp = c(rbind(from, to)), directed = TRUE)
    new <- eids == 0
    if (any(!new)) {
      graph <- igraph::set_edge_attr(
        graph,
        name = 'absPcorMin',
        index = eids[!new],
        value = abs(mapply(\(a, b) cor(Y[, a], Y[, b]), from[!new], to[!new]))
      )
    }
    # Add new edges
    from <- from[new]
    to <- to[new]
    if (length(from)) {
      abscors <- abs(mapply(\(a, b) cor(Y[, a], Y[, b]), from, to))
      graph <- igraph::add_edges(
        graph,
        edges = c(rbind(from, to)),
        attr = list(absPcorMin = abscors)
      )
    }
  }
  if (!allow_pgene_out) {
    npgenes <- setdiff(genes, c(pname, pgenes))
    if (length(npgenes)) {
      from <- rep(pgenes, each = length(npgenes))
      to <- rep(npgenes, times = length(pgenes))
      eids_to_delete <- setdiff(igraph::get_edge_ids(graph, vp = c(rbind(from, to)), directed = TRUE), 0)
      if (length(eids_to_delete)) graph <- igraph::delete_edges(graph, eids_to_delete)
    }
  }
  graph <- igraph::delete_vertices(graph, v = pname)
  return(graph)
}
