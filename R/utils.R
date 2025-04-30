#' Compute Pearson correlation matrix in parallel
#'
#' Efficiently computes the correlation matrix of a numeric matrix, using multiple cores for large matrices.
#'
#' @param mat Matrix.
#' @param ncores Number of cores.
#' @return Pearson correlation matrix.
#' @export
parallel_cor <- function(mat, ncores) {
  if ((nrow(mat) < 1e4) && (ncol(mat) < 100)) {
    return(cor(mat))
  }
  split_cols <- split(x = seq_len(ncol(mat)), f = cut(seq_len(ncol(mat)), ncores, labels = FALSE))
  cors <- parallel::mcmapply(function(cols) cor(mat[, cols], mat), split_cols, SIMPLIFY = F, mc.cores = ncores)
  return(do.call(rbind, cors))
}


#' Get Edge Index of Adjacency matrix
#'
#' Get edge index of adjacency matrix.
#'
#' @param G Adjacency matrix.
#' @return Edge index matrix.
#' @examples
#' G <- generate_complete_adj_matrix(genes = c('A', 'B'))
#' get_edge_index(G = G)
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
#' @param pMax P-value Matrix corresponding to edges in `G`.
#' @param chisqMin Chi-square statistics matrix corresponding to edges in `G`.
#' @return Directed igraph object with edge attributes `pMax` and `chisqMin`.
#' @export
adj2igraph <- function(G, pMax, chisqMin) {
  graph = igraph::graph_from_adjacency_matrix(G, mode = 'directed')
  # Delete unwanted edge attribute
  if ('weight' %in% igraph::edge_attr_names(graph)) graph = igraph::delete_edge_attr(graph, 'weight')
  # Add attributes to edges
  edges = do.call(rbind, strsplit(igraph::as_ids(igraph::E(graph)), '\\|'))
  igraph::E(graph)$pMax = pMax[edges]
  igraph::E(graph)$chisqMin = chisqMin[edges]
  return(graph)
}


#' Get separation set of two genes
#'
#' Extract separation set of two genes from graph skeleton
#'
#' @param skel See \code{\link{infer_skeleton_from_wildtype}} and \code{\link{infer_skeleton_from_perturbed}}.
#' @param g1 Gene name.
#' @param g2 Gene name.
#' @return Seperation set.
#' @export
get_sepset <- function(skel, g1, g2) {
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
