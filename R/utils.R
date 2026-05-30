#' Get Edge Index of Adjacency matrix
#'
#' Get edge index of adjacency matrix.
#'
#' @param G Adjacency matrix.
#' @return Edge index matrix.
#' @export
get_edge_index <- function(G) {
  G[lower.tri(G)] <- FALSE
  edge_index <- which(G, arr.ind = TRUE)
  rownames(edge_index) <- NULL
  return(edge_index)
}


#' Convert adjacency matrix to directed igraph with edge attributes
#'
#' Constructs a directed igraph object from an adjacency matrix, assigns `pMax` and `chisqMin` values to each edge.
#'
#' @param G Adjacency matrix.
#' @param pMax P-value matrix corresponding to edges in `G`.
#' @param chisqMin Chi-square statistics matrix corresponding to edges in `G`.
#' @param absPcorMin Partial correlation matrix corresponding to edges in `G`.
#' @param Threshold Threshold matrix corresponding to edges in `G`.
#' @param sampleSize sampleSize matrix corresponding to edges in `G`.
#' @return Directed igraph object with edge attributes `pMax` and `chisqMin`.
#' @export
adj2igraph <- function(G, pMax, chisqMin, absPcorMin, Threshold, sampleSize) {
  graph <- igraph::graph_from_adjacency_matrix(G, mode = 'directed')
  # Delete unwanted edge attribute
  if ('weight' %in% igraph::edge_attr_names(graph)) {
    graph <- igraph::delete_edge_attr(graph, 'weight')
  }
  # Add attributes to edges
  if (igraph::ecount(graph) > 0) {
    edges <- do.call(rbind, strsplit(igraph::as_ids(igraph::E(graph)), '\\|'))
    igraph::E(graph)$pMax <- pMax[edges]
    igraph::E(graph)$chisqMin <- chisqMin[edges]
    igraph::E(graph)$absPcorMin <- absPcorMin[edges]
    igraph::E(graph)$threshold <- Threshold[edges]
    igraph::E(graph)$n <- sampleSize[edges]
  }
  return(graph)
}


#' Get separation set of two genes
#'
#' Extract separation set of two genes from graph skeleton
#'
#' @param skel See \code{\link{infer_skeleton}}.
#' @param g1 Gene name.
#' @param g2 Gene name.
#' @return Separation set.
#' @export
get_sepset <- function(skel, g1, g2) {
  stopifnot(!is.null(skel$sepSet))
  nodes <- igraph::V(skel$graph)$name
  i <- match(g1, nodes)
  j <- match(g2, nodes)
  S <- skel$sepSet[[i]][[j]]
  if (is.null(S)) {
    return(NULL)
  } else {
    return(nodes[S])
  }
}
