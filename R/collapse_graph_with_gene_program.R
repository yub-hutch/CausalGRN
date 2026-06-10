#' Collapse graph with a gene program node
#'
#' Converts each incoming edge to the gene program node into incoming edges to
#' all member genes, then removes the gene program node.
#'
#' @param graph igraph object containing the gene program node.
#' @param pname Gene program node name.
#' @param pgenes Gene program member genes.
#' @param Y scRNA-seq normalized expression matrix with columns matching graph
#' vertex names.
#'
#' @return igraph object.
#' @export
collapse_graph_with_gene_program <- function(graph, pname, pgenes, Y) {
  stopifnot(inherits(graph, "igraph"))
  nodes <- igraph::V(graph)$name
  genes <- setdiff(nodes, pname)
  stopifnot(
    length(pname) == 1L,
    pname %in% nodes,
    is.character(pgenes),
    all(pgenes %in% genes),
    all(nodes %in% colnames(Y))
  )

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
        value = abs(mapply(
          \(a, b) stats::cor(Y[, a], Y[, b]),
          from[!new],
          to[!new]
        ))
      )
    }
    # Add new edges
    from <- from[new]
    to <- to[new]
    if (length(from)) {
      abscors <- abs(mapply(\(a, b) stats::cor(Y[, a], Y[, b]), from, to))
      graph <- igraph::add_edges(
        graph,
        edges = c(rbind(from, to)),
        attr = list(absPcorMin = abscors)
      )
    }
  }
  graph <- igraph::delete_vertices(graph, v = pname)
  return(graph)
}
