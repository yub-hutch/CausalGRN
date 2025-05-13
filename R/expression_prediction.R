#' Fit Expression Model
#'
#' Fit predictive models for gene expression using specified graph structure.
#'
#' @param Y Expression matrix (cells x genes).
#' @param group Vector indicating cell groups.
#' @param graph Graph structure ('mean', 'all', or igraph object).
#' @param ncores Number of cores for parallel computation.
#' @return List of fitted models, one per gene.
#' @export
fit_expression_model <- function(Y, group, graph, ncores) {
  p <- ncol(Y)
  genes <- colnames(Y)
  message(paste0('Fitting models for ', p, ' genes with ', nrow(Y), ' cells ...'))
  if (identical(graph, 'mean')) {
    graph <- igraph::graph_from_adjacency_matrix(matrix(FALSE, p, p, dimnames = list(genes, genes)), mode = 'directed')
  } else if (identical(graph, 'all')) {
    FC <- matrix(TRUE, p, p, dimnames = list(genes, genes))
    diag(FC) <- FALSE
    graph <- igraph::graph_from_adjacency_matrix(FC, mode = 'directed')
  }
  stopifnot('igraph' %in% class(graph))
  models <- pbmcapply::pbmclapply(genes, function(gene) {
    coef <- setNames(rep(0, p + 1), c('Intercept', genes))
    predictors <- igraph::neighbors(graph, gene, mode = 'in')$name
    samples <- which(group != gene)
    if (length(predictors) == 0) {
      coef['Intercept'] <- mean(Y[samples, gene])
    } else if (length(predictors) == 1) {
      coef[c('Intercept', predictors)] <- lm(Y[samples, gene] ~ Y[samples, predictors])$coefficients
    } else {
      cvfit <- glmnet::cv.glmnet(
        x = Y[samples, predictors],
        y = Y[samples, gene],
        family = "gaussian",
        nfolds = 5,
        nlambda = 100,
        alpha = 0,
        standardize = TRUE,
        intercept = TRUE,
        standardize.response = FALSE,
        parallel = FALSE
      )
      coef['Intercept'] <- cvfit$glmnet.fit$a0[cvfit$index["min", "Lambda"]]
      coef[predictors] <- cvfit$glmnet.fit$beta[, cvfit$index["min", "Lambda"]]
    }
    return(coef)
  }, mc.cores = ncores)
  return(models)
}

#' Apply Expression Model
#'
#' Predict gene expression for test cells using fitted models.
#'
#' @param models List of fitted models (See \code{\link{fit_expression_model}}).
#' @param Y Expression matrix (cells x genes).
#' @param group Vector indicating cell groups.
#' @param ncores Number of cores for parallel computation.
#' @return Predicted expression matrix (cells x genes).
#' @export
apply_expression_model <- function(models, Y, group, ncores) {
  stopifnot(identical(names(models), colnames(Y)))
  genes <- colnames(Y)
  pred_cell <- pbmcapply::pbmcmapply(function(gene, beta) {
    testX <- cbind(Intercept = 1, Y[, setdiff(genes, gene)])
    c(testX %*% beta[colnames(testX)])
  }, gene = names(models), beta = models, SIMPLIFY = TRUE, mc.cores = ncores)
  dimnames(pred_cell) <- dimnames(Y)
  kos <- sort(unique(group))
  pred <- matrix(NA, length(kos), length(genes), dimnames = list(kos, genes))
  for (ko in kos) {
    pred[ko, ] <- colMeans(pred_cell[group == ko, , drop = FALSE])
  }
  return(pred)
}


#' Apply Perturbation
#'
#' Simulate gene expression after perturbation using graph structure.
#'
#' @param ko_expressions Known KO gene expressions.
#' @param wt_expressions WT bulk gene expressions.
#' @param graph Graph structure (igraph object or 'all').
#' @param models List of fitted models (See \code{\link{fit_expression_model}}).
#' @param max_dist Maximum distance to propagate perturbation effect (default is Inf).
#' @return Predicted expression matrix (KO genes x all genes).
#' @export
apply_perturbation <- function(ko_expressions, wt_expressions, graph, models, max_dist = Inf) {
  impute_unknowns <- function(B, known_expressions) {
    stopifnot(all(diag(B[, -1]) == 0))
    stopifnot(identical(rownames(B), setdiff(colnames(B), 'Intercept')))
    stopifnot(all(names(known_expressions) %in% colnames(B)))
    genes <- rownames(B)
    U <- setdiff(genes, names(known_expressions))
    I_minus_B <- diag(length(U)) - B[U, U, drop = F]
    H <- 2 * t(I_minus_B) %*% I_minus_B
    eigenvals <- base::eigen(H, symmetric = TRUE)$values
    stopifnot(any(eigenvals > 0))
    stopifnot(max(eigenvals) / min(eigenvals) < 1e8)
    imputed <- base::solve(I_minus_B, B[U, names(known_expressions), drop = F] %*% known_expressions)[, 1]
    return(imputed)
  }
  p <- length(wt_expressions)
  genes <- names(wt_expressions)
  if (identical(graph, 'all')) {
    FC <- matrix(TRUE, p, p, dimnames = list(genes, genes))
    diag(FC) <- FALSE
    graph <- igraph::graph_from_adjacency_matrix(FC)
  }
  stopifnot('igraph' %in% class(graph))
  p <- length(wt_expressions)
  genes <- names(wt_expressions)
  ko_genes <- names(ko_expressions)
  stopifnot(setequal(igraph::V(graph)$name, genes))
  stopifnot(all(ko_genes %in% genes))
  B <- do.call(rbind, models) # Including intercept
  pred <- matrix(NA, length(ko_genes), p, dimnames = list(ko_genes, genes))
  for (ko_gene in ko_genes) {
    # Plug in provided expression of KO gene
    pred[ko_gene, ko_gene] <- ko_expressions[ko_gene]
    # For genes that are not descendants of KO gene, plug in WT bulk expression
    distances <- igraph::distances(graph, v = ko_gene, mode = 'out')[ko_gene, ]
    unchanged_genes <- names(which(distances >= max_dist + 1))
    if (length(unchanged_genes) > 0) {
      pred[ko_gene, unchanged_genes] <- wt_expressions[unchanged_genes]
    }
    # For descendants of KO gene, impute simultaneously since they can depend on each other
    known_expressions <- c(1, pred[ko_gene, c(ko_gene, unchanged_genes)])
    names(known_expressions) <- c('Intercept', c(ko_gene, unchanged_genes))
    if (length(known_expressions) < p + 1) {
      imputed <- impute_unknowns(B = B, known_expressions = known_expressions)
      pred[ko_gene, names(imputed)] <- imputed
    }
  }
  return(pred)
}
