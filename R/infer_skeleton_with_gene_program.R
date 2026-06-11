#' Infer graph skeleton with a gene program node
#'
#' Constructs a graph skeleton by performing iterative conditional independence
#' (CI) tests on all columns of the input matrices, including the gene program
#' node. Missing bidirectional skeleton edges are then added between the gene
#' program node and its member genes.
#'
#' @param count scRNA-seq count matrix (cells x nodes, including gene program).
#' @param Y scRNA-seq normalized expression matrix (cells × nodes, including
#' gene program). For example, log1p(total UMI corrected count).
#' @param alpha Significance level for CI tests.
#' @param pname Name of the gene program node, such as \code{"PC1"}.
#' @param pgenes Character vector of gene program member genes.
#' @param min_abspcor Minimum absolute value of partial correlation for kept edges.
#' @param ncores Number of CPU cores for parallel processing.
#' @param G Optional initial adjacency matrix (defaults to a fully connected
#' graph without self-loops).
#' @param max_order Maximum conditioning set size (0 or 1, default is 1).
#' @param max_thr Maximum threshold for conditional variable (default is 10).
#' @param min_n1 Minimum number of samples satisfying Yk > selected threshold
#' (default is 1000).
#' @param min_n2 Minimum number of samples satisfying Yk > selected threshold,
#' Yi > 0, and Yj > 0 (default is 200).
#'
#' @return An igraph object representing the inferred skeleton.
#' @export
infer_skeleton_with_gene_program <- function(
    count, Y, alpha, min_abspcor, pname, pgenes, ncores, G = NULL,
    max_order = 1, max_thr = 10, min_n1 = 1000, min_n2 = 200
) {
  nodes <- colnames(Y)
  .check_gene_program_nodes(pname = pname, pgenes = pgenes, nodes = nodes)

  graph <- infer_skeleton(
    count = count,
    Y = Y,
    alpha = alpha,
    min_abspcor = min_abspcor,
    ncores = ncores,
    G = G,
    max_order = max_order,
    max_thr = max_thr,
    min_n1 = min_n1,
    min_n2 = min_n2,
    sepset = FALSE
  )$graph
  # Add edges connecting gene program and program genes.
  pgenes_to_add <- setdiff(pgenes, igraph::neighbors(graph, v = pname)$name)
  if (length(pgenes_to_add)) {
    from <- c(rep(pname, length(pgenes_to_add)), pgenes_to_add)
    to <- c(pgenes_to_add, rep(pname, length(pgenes_to_add)))
    abscors <- abs(c(stats::cor(Y[, pname], Y[, pgenes_to_add])))
    graph <- igraph::add_edges(
      graph,
      edges = c(rbind(from, to)),
      attr = list(absPcorMin = c(abscors, abscors))
    )
  }
  return(graph)
}
