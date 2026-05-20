#' Run PC Algorithm
#'
#' Runs PC algorithm and return an igraph object.
#'
#' @param wt scRNA-seq matrix of wild-type cells.
#' @param alpha Significance level for conditional independence tests.
#' @param ncores Number of cores to use for parallel computation (default is 1).
#' @return igraph object.
#' @export
run_pc <- function(wt, alpha, ncores = 1) {
  .check_expression_matrix(wt, arg = "wt")
  if (
    length(alpha) != 1L || !is.numeric(alpha) ||
      is.na(alpha) || !is.finite(alpha) || alpha <= 0 || alpha >= 1
  ) {
    stop("'alpha' must be a single finite number between 0 and 1.", call. = FALSE)
  }
  .check_ncores(ncores)

  # Fit
  raw <- pcalg::pc(
    suffStat = list(C = cor(wt), n = nrow(wt)),
    indepTest = pcalg::gaussCItest,
    alpha = alpha,
    labels = colnames(wt),
    u2pd = 'relaxed',
    skel.method = 'stable.fast',
    conservative = FALSE,
    maj.rule = FALSE,
    solve.confl = TRUE,
    numCores = ncores,
    verbose = FALSE
  )

  # With stable.fast, the lower triangle and diagonal of pMax are filled with -1.
  diag(raw@pMax) <- 1
  raw@pMax[lower.tri(raw@pMax)] <- t(raw@pMax)[lower.tri(raw@pMax)]
  rownames(raw@pMax) <- colnames(raw@pMax) <- colnames(wt)

  # Get igraph object
  graph <- igraph::graph_from_graphnel(raw@graph)
  if (!igraph::is_directed(graph)) {
    graph <- igraph::as_directed(graph, mode = 'mutual')
  }

  # Delete unwanted edge attribute
  if ('weight' %in% igraph::edge_attr_names(graph)) {
    graph <- igraph::delete_edge_attr(graph, name = 'weight')
  }

  # Add pMax and score attributes to edges.
  if (igraph::ecount(graph) == 0L) {
    return(graph)
  }
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
  .check_expression_matrix(wt, arg = "wt")

  # Fit
  raw <- pcalg::ges(
    score = methods::new('GaussL0penObsScore', wt),
    verbose = verbose
  )

  # Convert to igraph object
  amat <- as(as(raw$essgraph, 'graphNEL'), 'matrix') > 0
  graph <- igraph::graph_from_adjacency_matrix(
    amat,
    mode = 'directed',
    weighted = NULL,
    diag = FALSE
  )
  igraph::E(graph)$score <- 1
  return(graph)
}


#' Run GIES Algorithm
#'
#' Run GIES algorithm and return an igraph object.
#'
#' @param wt scRNA-seq matrix of wild-type cells.
#' @param pts List of scRNA-seq data matrices of perturbed cells, with names
#'   representing perturbed gene.
#' @param verbose Whether to print training details (default is FALSE).
#' @return igraph object.
#' @export
run_gies <- function(wt, pts, verbose = FALSE) {
  .check_expression_matrix(wt, arg = "wt")
  if (!is.list(pts) || length(pts) == 0L) {
    stop(
      "'pts' must be a non-empty list of perturbed expression matrices.",
      call. = FALSE
    )
  }
  if (is.null(names(pts)) || anyNA(names(pts)) || any(names(pts) == "")) {
    stop("'pts' must be named by perturbed gene.", call. = FALSE)
  }
  if (anyDuplicated(names(pts))) {
    stop("'pts' names must be unique.", call. = FALSE)
  }

  # Prepare inputs
  genes <- colnames(wt)
  targets <- names(pts)
  if (!all(targets %in% genes)) {
    stop("'pts' names must all be columns of 'wt'.", call. = FALSE)
  }
  for (i in seq_along(pts)) {
    pt <- pts[[i]]
    .check_expression_matrix(pt, arg = paste0("pts[[", i, "]]"))
    if (!identical(colnames(pt), genes)) {
      stop("Each matrix in 'pts' must have the same columns as 'wt'.", call. = FALSE)
    }
  }
  combined <- rbind(wt, Reduce(rbind, pts))

  # Specify which genes are intervened in each configuration
  targets_gies <- c(list(integer(0)), as.list(match(targets, genes)))

  # Specify which configuration each cell belongs to
  target_index_gies <- c(
    rep(1, nrow(wt)),
    unlist(lapply(seq_along(pts), \(i) rep(i + 1, nrow(pts[[i]]))))
  )
  score <- methods::new(
    'GaussL0penIntScore',
    data = combined,
    targets = targets_gies,
    target.index = target_index_gies,
    nodes = genes
  )

  # Fit
  raw <- pcalg::gies(score, verbose = verbose)

  # Convert to igraph object
  amat <- as(as(raw$essgraph, 'graphNEL'), 'matrix') > 0
  graph <- igraph::graph_from_adjacency_matrix(
    amat,
    mode = 'directed',
    weighted = NULL,
    diag = FALSE
  )
  igraph::E(graph)$score <- 1
  return(graph)
}


#' Run GENIE3 Algorithm
#'
#' Run GENIE3 algorithm and return an igraph object.
#'
#' @param wt scRNA-seq matrix of wild-type cells.
#' @param ncores Number of cores.
#' @param verbose Whether to print training details (default is FALSE).
#' @return igraph object.
#' @export
run_genie3 <- function(wt, ncores, verbose = FALSE) {
  if (!requireNamespace("GENIE3", quietly = TRUE)) {
    stop("Package 'GENIE3' must be installed to use run_genie3().", call. = FALSE)
  }
  .check_expression_matrix(wt, arg = "wt")
  .check_ncores(ncores)

  w <- GENIE3::GENIE3(t(wt), nCores = ncores, verbose = verbose)
  w <- w[colnames(wt), colnames(wt)]
  graph <- igraph::graph_from_adjacency_matrix(
    abs(w),
    mode = 'directed',
    weighted = TRUE,
    diag = FALSE
  )
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
  .check_expression_matrix(Y)

  w <- cor(Y)
  diag(w) <- 0
  graph <- igraph::graph_from_adjacency_matrix(
    abs(w),
    mode = 'directed',
    weighted = TRUE,
    diag = FALSE
  )
  igraph::E(graph)$score <- igraph::E(graph)$weight
  graph <- igraph::delete_edge_attr(graph, name = 'weight')
  return(graph)
}


#' Format GRNBoost2 Output
#'
#' Reads GRNBoost2 output csv file and returns an igraph object.
#'
#' @param file GRNBoost2 output.
#' @param genes Gene names.
#' @return igraph object.
#' @export
format_grnboost2 <- function(file, genes) {
  df <- read.csv(file, header = TRUE, stringsAsFactors = FALSE)
  required_cols <- c("TF", "target", "importance")
  if (!all(required_cols %in% colnames(df))) {
    stop("'file' must contain columns: TF, target, importance.", call. = FALSE)
  }
  if (
    !is.character(genes) || length(genes) == 0L ||
      anyNA(genes) || any(genes == "") || anyDuplicated(genes)
  ) {
    stop("'genes' must be a non-missing unique character vector.", call. = FALSE)
  }
  if (!is.numeric(df$importance) || any(!is.finite(df$importance))) {
    stop("'importance' must contain only finite numeric values.", call. = FALSE)
  }
  if (!all(df$TF %in% genes) || !all(df$target %in% genes)) {
    stop(
      "All TF and target values in 'file' must be included in 'genes'.",
      call. = FALSE
    )
  }

  raw <- with(
    df,
    tapply(importance, INDEX = list(TF, target), FUN = mean, default = 0)
  )
  w <- matrix(0, length(genes), length(genes), dimnames = list(genes, genes))
  w[rownames(raw), colnames(raw)] <- raw
  graph <- igraph::graph_from_adjacency_matrix(
    abs(w),
    mode = 'directed',
    weighted = TRUE,
    diag = FALSE
  )
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
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("Package 'glmnet' must be installed to use run_lasso().", call. = FALSE)
  }
  .check_expression_matrix(Y)
  .check_ncores(ncores)
  if (
    length(nfold) != 1L || !is.numeric(nfold) ||
      is.na(nfold) || !is.finite(nfold) ||
      nfold < 2 || nfold != round(nfold)
  ) {
    stop("'nfold' must be a single integer >= 2.", call. = FALSE)
  }
  if (ncol(Y) <= 2L) {
    stop("'Y' must contain more than two genes.", call. = FALSE)
  }

  genes <- colnames(Y)
  coef_list <- .parallel_lapply(
    genes,
    function(g) {
      coef <- setNames(rep(NA, ncol(Y)), colnames(Y))
      coef[g] <- 0
      x <- Y[, setdiff(colnames(Y), g), drop = FALSE]
      y <- Y[, g]
      if (ncol(x) == 1) {
        coef[colnames(x)] <- stats::lm(y ~ x)$coefficients['x']
      } else {
        coef[colnames(x)] <- tryCatch({
          cvfit <- glmnet::cv.glmnet(
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
      coef
    },
    ncores = ncores,
    preschedule = FALSE,
    export = c("Y", "nfold")
  )
  w <- do.call(cbind, coef_list)
  colnames(w) <- genes
  graph <- igraph::graph_from_adjacency_matrix(
    abs(w),
    mode = 'directed',
    weighted = TRUE,
    diag = FALSE
  )
  igraph::E(graph)$score <- igraph::E(graph)$weight
  graph <- igraph::delete_edge_attr(graph, name = 'weight')
  return(graph)
}
