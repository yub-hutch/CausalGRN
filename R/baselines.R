#' Run PC Algorithm
#'
#' Runs PC algorithm and return an igraph object.
#'
#' @param wt scRNA-seq matrix of wild-type cells.
#' @param alpha Significance level for conditional independence tests.
#' @return igraph object.
#' @export
run_pc <- function(wt, alpha) {
  # Fit
  raw <- pcalg::pc(
    suffStat = list(C = cor(wt), n = nrow(wt)),
    indepTest = pcalg::gaussCItest,
    alpha = alpha,
    labels = colnames(wt),
    u2pd = 'relaxed',
    skel.method = 'stable.fast',
    conservative = F,
    maj.rule = F,
    solve.confl = T,
    numCores = 1,
    verbose = F
  )
  # Format pMax
  # When skel.method is stable.fast, the lower triangle & Diagonal of pMax are filled with -1
  diag(raw@pMax) <- 1
  raw@pMax[lower.tri(raw@pMax)] <- t(raw@pMax)[lower.tri(raw@pMax)]
  rownames(raw@pMax) <- colnames(raw@pMax) <- colnames(wt)
  # Get igraph object
  graph <- igraph::graph_from_graphnel(raw@graph)
  if (!igraph::is_directed(graph)) {
    graph <- igraph::as.directed(graph, mode = 'mutual')
  }
  # Delete unwanted edge attribute
  if ('weight' %in% igraph::edge_attr_names(graph)) {
    graph <- igraph::delete_edge_attr(graph, name = 'weight')
  }
  # Add pMax & score attribute to edges
  edges <- do.call(rbind, strsplit(igraph::as_ids(igraph::E(graph)), '\\|'))
  igraph::E(graph)$pMax <- raw@pMax[edges]
  igraph::E(graph)$score <- 1 - raw@pMax[edges]
  return(graph)
}


#' Run GES Algorithm
#'
#' Run a GES algorithm and return an igraph object.
#'
#' @param wt scRNA-seq matrix of wild-type cells.
#' @param verbose Whether to print training details (default is FALSE).
#' @return igraph object.
#' @export
run_ges <- function(wt, verbose = FALSE) {
  # Fit
  raw <- pcalg::ges(score = new('GaussL0penObsScore', wt), verbose = verbose)
  # Convert to igraph object
  amat <- as(as(raw$essgraph, 'graphNEL'), 'matrix') > 0
  graph <- igraph::graph_from_adjacency_matrix(amat, mode = 'directed', weighted = NULL, diag = FALSE)
  igraph::E(graph)$score <- 1
  return(graph)
}


#' Run GIES Algorithm
#'
#' Run GIES algorithm and return an igraph object.
#'
#' @param wt scRNA-seq matrix of wild-type cells.
#' @param pts List of scRNA-seq data matrices of perturbed cells, with names representing perturbed gene.
#' @param verbose Whether to print training details (default is FALSE).
#' @return igraph object.
#' @export
run_gies <- function(wt, pts, verbose = FALSE) {
  # Prepare inputs
  genes <- colnames(wt)
  targets <- names(pts)
  combined <- rbind(wt, Reduce(rbind, pts))
  # Specify which genes are intervened in each configuration
  targets_gies <- c(list(integer(0)), as.list(match(targets, genes)))
  # Specify which configuration each cell belongs to
  target_index_gies <- c(rep(1, nrow(wt)), unlist(lapply(seq_along(pts), function(i) rep(i + 1, nrow(pts[[i]])))))
  score <- new('GaussL0penIntScore', data = combined, targets = targets_gies, target.index = target_index_gies, nodes = genes)
  # Fit
  raw <- pcalg::gies(score, verbose = verbose)
  # Convert to igraph object
  amat <- as(as(raw$essgraph, 'graphNEL'), 'matrix') > 0
  graph <- igraph::graph_from_adjacency_matrix(amat, mode = 'directed', weighted = NULL, diag = FALSE)
  igraph::E(graph)$score <- 1
  return(graph)
}


#' Run GENIE3 Algorithm
#'
#' Run GENIE3 algorithm and return an importance matrix.
#'
#' @param wt scRNA-seq matrix of wild-type cells.
#' @param ncores Number of cores.
#' @param verbose Whether to print training details (default is FALSE).
#' @return igraph object.
#' @export
run_genie3 <- function(wt, ncores, verbose = FALSE) {
  w <- GENIE3::GENIE3(t(wt), nCores = ncores, verbose = verbose)
  w <- w[colnames(wt), colnames(wt)]
  graph <- igraph::graph_from_adjacency_matrix(abs(w), mode = 'directed', weighted = TRUE, diag = FALSE)
  igraph::E(graph)$score <- igraph::E(graph)$weight
  graph <- igraph::delete_edge_attr(graph, name = 'weight')
  return(graph)
}


#' Build correlation based GRN
#'
#' Build GRN based on marginal correlations between genes.
#'
#' @param Y Matrix of normalized scRNA-seq data.
#' @return igraph object.
#' @export
run_cor <- function(Y) {
  w <- cor(Y)
  diag(w) <- 0
  graph <- igraph::graph_from_adjacency_matrix(abs(w), mode = 'directed', weighted = TRUE, diag = FALSE)
  igraph::E(graph)$score <- igraph::E(graph)$weight
  graph <- igraph::delete_edge_attr(graph, name = 'weight')
  return(graph)
}


#' Format GRNBoost2 Output
#'
#' Reads GRNBoost2 output csv file and returns importance matrix.
#'
#' @param file GRNBoost2 output.
#' @param genes Gene names.
#' @return igraph object.
#' @export
format_grnboost2 <- function(file, genes) {
  df <- read.csv(file, header = TRUE, stringsAsFactors = FALSE)
  raw <- with(df, tapply(importance, INDEX = list(TF, target), FUN = mean, default = 0))
  w <- matrix(0, length(genes), length(genes), dimnames = list(genes, genes))
  w[rownames(raw), colnames(raw)] <- raw
  graph <- igraph::graph_from_adjacency_matrix(abs(w), mode = 'directed', weighted = TRUE, diag = FALSE)
  igraph::E(graph)$score <- igraph::E(graph)$weight
  graph <- igraph::delete_edge_attr(graph, name = 'weight')
  return(graph)
}


#' Infer GRN with Lasso
#'
#' Infers GRN with lasso, treating all cells as homogeneous.
#'
#' @param Y Matrix of normalized scRNA-seq data.
#' @param ncores Number of cores to use for parallel computation.
#' @param nfold Number of folds for cross-validation (default is 5).
#' @return igraph object.
#' @export
run_lasso <- function(Y, ncores, nfold = 5) {
  stopifnot(ncol(Y) > 2)
  w <- pbmcapply::pbmcmapply(function(g) {
    coef <- setNames(rep(NA, ncol(Y)), colnames(Y))
    coef[g] <- 0
    x <- Y[, setdiff(colnames(Y), g), drop = FALSE]
    y <- Y[, g]
    if (ncol(x) == 1) {
      coef[colnames(x)] <- lm(y ~ x)$coefficients['x']
    } else {
      coef[colnames(x)] <- tryCatch({
        cvfit = glmnet::cv.glmnet(
          x = x,
          y = y,
          family = 'gaussian',
          nfolds = nfold,
          nlambda = 100,
          alpha = 1,
          standardize = TRUE,
          intercept = TRUE,
          standardize.response = FALSE,
          parallel = FALSE
        )
        cvfit$glmnet.fit$beta[, cvfit$index['min', 'Lambda']]
      }, error = function(e) rep(0, ncol(x)))
    }
    return(coef)
  }, colnames(Y), SIMPLIFY = TRUE, mc.cores = ncores, mc.preschedule = FALSE)
  graph <- igraph::graph_from_adjacency_matrix(abs(w), mode = 'directed', weighted = TRUE, diag = FALSE)
  igraph::E(graph)$score <- igraph::E(graph)$weight
  graph <- igraph::delete_edge_attr(graph, name = 'weight')
  return(graph)
}
