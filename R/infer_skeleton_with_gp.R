#' Infer graph skeleton with gene program
#'
#' Constructs graph skeleton by performing iterative conditional independence (CI) tests on gene expression data.
#' Edges are removed based on statistical independence.
#'
#' @param count scRNA-seq count matrix (cells x genes).
#' @param Y scRNA-seq normalized expression matrix (cells × genes, including gene program). For example, log1p(total UMI corrected count).
#' @param alpha Significance level for CI tests.
#' @param pname Gene program name in \code{Y}.
#' @param pgenes Gene program genes in \code{Y}.
#' @param min_abspcor Minimum absolute value of partial correlation for kept edges.
#' @param ncores Number of CPU cores for parallel processing.
#' @param G Optional initial adjacency matrix (defaults to a fully connected graph without self-loops).
#' @param max_order Maximum conditioning set size (0 or 1, default is 1).
#' @param max_thr Maximum threshold for conditional variable (default is 10).
#' @param min_n1 Minimum number of samples satisfying Yk > selected threshold (default is 2000).
#' @param min_n2 Minimum number of samples satisfying Yk > selected threshold, Yi > 0, and Yj > 0 (default is 400).
#' @param max_nchildren Maximum number of children a node can have (default is Inf).
#' @param max_nparent Maximum number of parents a node can have (default is Inf).
#' @param sepset Return separation set or not (default is \code{TRUE}).
#'
#' @return A list with:
#' \describe{
#'   \item{graph}{An igraph object representing the inferred skeleton with edge attributes `pMax`, `chisqMin`, and `absPcorMin`.}
#'   \item{sepSet}{List of separation sets for each node pair.}
#' }
#' @export
infer_skeleton_with_gp <- function(
    count, Y, alpha, min_abspcor, pname, pgenes, ncores, G = NULL, max_order = 1, max_thr = 10, min_n1 = 2000, min_n2 = 400,
    max_nchildren = Inf, max_nparent = Inf, sepset = TRUE
) {
  genes = colnames(Y)
  stopifnot(
    pname %in% genes,
    all(pgenes %in% genes)
  )
  res <- infer_skeleton(
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
    max_nchildren = max_nchildren,
    max_nparent = max_nparent,
    sepset = sepset
  )
  # Add edges connecting gene program and program genes
  pgenes_to_add <- setdiff(pgenes, igraph::neighbors(res$graph, v = pname)$name)
  if (length(pgenes_to_add)) {
    edges_to_add <- as.vector(rbind(pname, pgenes_to_add))
    res$graph <- igraph::add_edges(res$graph, edges = edges_to_add)
    if (sepset) {
      p_idx <- match(pname, genes)
      pgene_indices_to_update <- match(pgenes_to_add, genes)
      for (pgene_idx in pgene_indices_to_update) {
        res$sepSet[[p_idx]][[pgene_idx]] <- res$sepSet[[pgene_idx]][[p_idx]] <- NULL
      }
    }
  }
  return(res)
}
