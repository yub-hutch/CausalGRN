#' Predict KO Effect Using the Mean Perturbation Effect Baseline
#'
#' This function calculates the average perturbation effect (delta from wild-type)
#' from a training set and uses it to predict the effect for a test set.
#'
#' @param Y A numeric matrix of expression data (cells x genes).
#' @param group A character or factor vector indicating the group for each cell in Y.
#' @param train_kos A character vector of the training knockout group names.
#' @param test_kos A character vector of the test knockout group names.
#' @param wt_name A character string specifying the wild-type group name (default is "WT").
#'
#' @return A numeric matrix of predicted delta values (rows are test KOs, columns are genes).
#' @export
predict_mean_effect <- function(Y, group, train_kos, test_kos, wt_name = "WT") {

  # 1. Perform aggressive data checks
  stopifnot(wt_name %in% group)
  stopifnot(all(train_kos %in% group))
  stopifnot(all(test_kos %in% group))

  # 2. Calculate the mean wild-type expression vector
  wt_expressions <- colMeans(Y[group == wt_name, , drop = FALSE])

  # 3. Calculate the "Mean Delta" from the training data
  train_delta_matrix <- vapply(train_kos, function(ko) {
    colMeans(Y[group == ko, , drop = FALSE]) - wt_expressions
  }, FUN.VALUE = numeric(ncol(Y)))

  # Average all the delta vectors together
  mean_delta <- rowMeans(train_delta_matrix)

  # 4. Create the final output matrix of predicted deltas
  # Each row is the identical mean delta vector.
  pred_matrix <- t(replicate(length(test_kos), mean_delta))
  rownames(pred_matrix) <- test_kos

  return(pred_matrix)
}


#' Fit Gene Expression Model
#'
#' Fits a predictive model for each gene's expression based on a given graph
#' structure and regression method.
#'
#' @param Y A numeric matrix of expression data (cells x genes).
#' @param group A character or factor vector indicating the group for each cell in Y.
#' @param graph An igraph object or the character string 'all' specifying the
#'   regulatory structure.
#' @param ncores The number of cores for parallel computation.
#' @param method The regression method to use. One of 'lm', 'lasso', or 'ridge'.
#' @return A numeric matrix `B` where rows are sources ('Intercept' and genes)
#'   and columns are targets (genes).
#' @export
fit_expression_model <- function(Y, group, graph, ncores, method = c('lm', 'lasso', 'ridge')) {
  # --- 1. Parameter and Graph Setup ---
  p <- ncol(Y)
  genes <- colnames(Y)
  method <- match.arg(method)

  # --- UPDATED: More robust parameter checking ---
  if (inherits(graph, "igraph")) {
    stopifnot(setequal(genes, igraph::V(graph)$name))
    stopifnot(!igraph::any_loop(graph))
  } else {
    stopifnot(identical(graph, "all"))
  }

  message(paste0('Fitting models for ', p, ' genes with ', nrow(Y), ' cells using method: ', method, ' ...'))

  if (identical(graph, 'all')) {
    adj_matrix <- matrix(TRUE, p, p, dimnames = list(genes, genes))
    diag(adj_matrix) <- FALSE
    graph <- igraph::graph_from_adjacency_matrix(adj_matrix, mode = 'directed')
  }

  # --- 2. Fit Models in Parallel ---
  model_list <- pbmcapply::pbmclapply(genes, function(gene) {
    # Initialize the full coefficient vector for this target gene
    coef_vector <- setNames(rep(0, p + 1), c('Intercept', genes))

    predictors <- igraph::neighbors(graph, gene, mode = 'in')$name
    samples <- which(group != gene)

    if (length(predictors) == 0) {
      coef_vector['Intercept'] <- mean(Y[samples, gene])
    } else {
      # --- BUG FIX: Use explicit assignment with `predictors` to avoid name sanitization issues ---
      if (length(predictors) == 1 || method == 'lm') {
        fit <- lm(Y[samples, gene] ~ Y[samples, predictors])

        if (!any(is.na(stats::coef(fit)))) {
          # Use the known, correct names for assignment
          coef_vector[c('Intercept', predictors)] <- stats::coef(fit)
        } else {
          warning(paste("lm failed for gene", gene, "- fitting intercept only."))
          coef_vector['Intercept'] <- mean(Y[samples, gene])
        }

      } else { # Use glmnet for lasso or ridge
        alpha_val <- ifelse(method == 'lasso', 1, 0)
        cvfit <- glmnet::cv.glmnet(
          x = Y[samples, predictors, drop = FALSE],
          y = Y[samples, gene],
          family = "gaussian",
          nfolds = 5,
          alpha = alpha_val,
          intercept = TRUE
        )
        fit_coefs <- as.matrix(stats::coef(cvfit, s = "lambda.min"))

        # Use the known, correct names for assignment
        coef_vector['Intercept'] <- fit_coefs[1, 1]
        coef_vector[predictors] <- fit_coefs[-1, 1]
      }
    }

    return(coef_vector)
  }, mc.cores = ncores)

  # --- 3. Assemble and Return the Final B Matrix ---
  B <- do.call(cbind, model_list)
  colnames(B) <- genes

  stopifnot(!any(is.na(B)))
  return(B)
}


#' Predict Oracle Perturbation Effect (Delta)
#'
#' Predicts the perturbation effect (delta from wild-type) for test cells using
#' a fitted B matrix, where the true expression values of the predictors are
#' known (oracle setting).
#'
#' @param B A numeric matrix of regulatory coefficients where rows are sources
#'   ('Intercept' and genes) and columns are targets (genes). This is the output
#'   of \code{\link{fit_expression_model}}.
#' @param test_Y A numeric matrix of expression data (cells x genes) for which to
#'   make predictions.
#' @param test_group A character or factor vector indicating the group for each
#'   cell in `test_Y`.
#' @param wt_expressions A named numeric vector of the wild-type expression levels
#'   for each gene.
#' @return A numeric matrix of predicted delta values (rows are unique KOs,
#'   columns are genes).
#' @export
predict_oracle_effect <- function(B, test_Y, test_group, wt_expressions) {
  # --- 1. Input Validation ---
  genes <- colnames(test_Y)
  stopifnot(
    is.matrix(B), is.matrix(test_Y),
    is.numeric(wt_expressions), !is.null(names(wt_expressions)),
    identical(colnames(B), genes),
    identical(rownames(B), c('Intercept', genes)),
    identical(names(wt_expressions), genes)
  )

  # --- 2. Perform Oracle Prediction for Absolute Expression ---
  # Prepare the full input matrix by adding an intercept column to the test data.
  X_full <- cbind(Intercept = 1, test_Y)

  # The core prediction is a single, efficient matrix multiplication.
  pred_cell_abs <- X_full %*% B

  # --- 3. Aggregate Cell-Level Predictions to KO-Level ---
  kos <- sort(unique(as.character(test_group)))

  # Aggregate by taking the mean absolute prediction for all cells in a KO group.
  pred_ko_abs <- vapply(kos, function(ko) {
    colMeans(pred_cell_abs[test_group == ko, , drop = FALSE])
  }, FUN.VALUE = numeric(length(genes)))

  pred_ko_abs <- t(pred_ko_abs)

  # --- 4. Convert Absolute Predictions to Delta Predictions ---
  # Subtract the WT expression from each row to get the predicted effect.
  pred_delta <- pred_ko_abs - matrix(wt_expressions, nrow = nrow(pred_ko_abs), ncol = ncol(pred_ko_abs), byrow = TRUE)

  return(pred_delta)
}


#' Impute Unknown Perturbation Effects (Delta)
#'
#' A helper function that solves the linear system Δx = BΔx to find the
#' steady-state deltas for a set of unknown genes, given a set of known deltas.
#'
#' @param B_propagator A square, gene-by-gene regulatory matrix (targets x sources).
#' @param known_deltas A named numeric vector of the initial, known deltas.
#' @return A named numeric vector of the imputed delta values for the unknown genes.
#' @noRd
impute_effect_deltas <- function(B_propagator, known_deltas) {
  genes <- rownames(B_propagator)
  unknown_genes <- setdiff(genes, names(known_deltas))

  if (length(unknown_genes) == 0) {
    return(numeric(0))
  }

  B_UU <- B_propagator[unknown_genes, unknown_genes, drop = FALSE]

  # Stricter Stability Check with Tolerance
  epsilon <- 1e-8
  eigen_vals <- base::eigen(B_UU, only.values = TRUE)$values
  if (any(abs(eigen_vals) >= 1 - epsilon)) {
    warning("System is unstable or nearly unstable (max |eigenvalue| >= 1 - epsilon). Results may be unreliable.")
  }

  # Calculate the driving force from the known genes ONLY. The intercept is not part of the delta.
  driving_force_known <- B_propagator[unknown_genes, names(known_deltas), drop = FALSE] %*% known_deltas

  # Solve for the deltas of the unknown genes
  I_minus_B_UU <- diag(length(unknown_genes)) - B_UU
  imputed_deltas <- base::solve(I_minus_B_UU, driving_force_known)

  # Ensure the output is a simple named vector
  return(imputed_deltas[,1])
}


#' Predict Standard Perturbation Effect (Delta)
#'
#' Predicts the perturbation effect (delta from wild-type) for a set of KOs
#' by propagating the initial perturbation through the network structure.
#'
#' @param B A numeric matrix of regulatory coefficients where rows are sources
#'   ('Intercept' and genes) and columns are targets (genes). This is the output
#'   of \code{\link{fit_expression_model}}.
#' @param ko_expressions A named numeric vector of the expression levels for the
#'   perturbed genes (e.g., 0 for a full knockout).
#' @param wt_expressions A named numeric vector of the wild-type expression levels
#'   for each gene.
#' @param graph An optional igraph object. If NULL (default), a functional graph is
#'   derived from B. If provided, this "roadmap" graph is used for distance
#'   calculations.
#' @param max_dist The maximum distance from the KO gene for an effect to be
#'   propagated. Default is `Inf` (full propagation).
#' @return A numeric matrix of predicted delta values (rows are KOs, columns are genes).
#' @export
predict_standard_effect <- function(B, ko_expressions, wt_expressions, graph = NULL, max_dist = Inf) {
  # --- 1. Input Validation and Setup ---
  genes <- names(wt_expressions)
  ko_genes <- names(ko_expressions)
  stopifnot(
    is.matrix(B), is.numeric(wt_expressions), is.numeric(ko_expressions),
    !is.null(names(wt_expressions)), !is.null(names(ko_expressions)),
    identical(colnames(B), genes),
    identical(rownames(B), c('Intercept', genes)),
    identical(names(wt_expressions), genes),
    all(ko_genes %in% genes)
  )

  # --- 2. Prepare B Matrix and Determine Propagation Graph ---
  B_propagator <- t(B[genes, , drop = FALSE])

  if (is.null(graph)) {
    # Default behavior: create functional graph from B matrix
    adj_matrix_functional <- t(B_propagator != 0)
    graph_for_distances <- igraph::graph_from_adjacency_matrix(adj_matrix_functional, mode = 'directed')
  } else {
    # Use the provided "roadmap" graph
    stopifnot('igraph' %in% class(graph), setequal(igraph::V(graph)$name, genes))
    graph_for_distances <- graph

    # --- BUG FIX: Correctly check that the functional graph is a subset of the roadmap graph ---
    adj_matrix_functional <- t(B_propagator != 0)
    graph_functional <- igraph::graph_from_adjacency_matrix(adj_matrix_functional, mode = 'directed')

    # The number of edges in the difference between the two graphs must be 0.
    stopifnot(
      igraph::ecount(igraph::difference(graph_functional, graph_for_distances)) == 0,
      "Functional graph from B contains edges not present in the provided roadmap graph."
    )
  }

  # --- 3. Main Prediction Loop ---
  pred_delta_list <- lapply(ko_genes, function(ko_gene) {

    # --- a. Define the Initial State (Known Deltas) ---
    known_deltas <- list()

    known_deltas[[ko_gene]] <- ko_expressions[ko_gene] - wt_expressions[ko_gene]

    # Identify non-descendants using the chosen graph for distances.
    distances <- igraph::distances(graph_for_distances, v = ko_gene, mode = 'out')[1, ]
    unchanged_genes <- names(which(is.infinite(distances) | distances > max_dist))
    for (ug in unchanged_genes) {
      known_deltas[[ug]] <- 0
    }

    clean_names <- names(known_deltas)
    known_deltas <- unlist(known_deltas, use.names = FALSE)
    names(known_deltas) <- clean_names

    # --- b. Solve for the Unknown Gene Deltas ---
    imputed_deltas <- impute_effect_deltas(
      B_propagator = B_propagator,
      known_deltas = known_deltas
    )

    # --- c. Assemble the Full Delta Vector for this KO ---
    full_delta_vector <- setNames(rep(0, length(genes)), genes)
    full_delta_vector[names(known_deltas)] <- known_deltas
    if (length(imputed_deltas) > 0) {
      full_delta_vector[names(imputed_deltas)] <- imputed_deltas
    }

    return(full_delta_vector)
  })

  # --- 4. Assemble and Return the Final Delta Matrix ---
  pred_matrix <- do.call(rbind, pred_delta_list)
  rownames(pred_matrix) <- ko_genes

  return(pred_matrix)
}
