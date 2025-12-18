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
.impute_deltas <- function(B_propagator, known_deltas) {
  genes <- rownames(B_propagator)
  unknown_genes <- setdiff(genes, names(known_deltas))
  if (length(unknown_genes) == 0) return(numeric(0))

  B_UU <- B_propagator[unknown_genes, unknown_genes, drop = FALSE]

  epsilon <- 1e-8
  eigen_vals <- base::eigen(B_UU, only.values = TRUE)$values
  if (any(abs(eigen_vals) >= 1 - epsilon)) {
    warning("System may be unstable (max |eigenvalue| >= 1 - epsilon).")
  }

  driving_force <- B_propagator[unknown_genes, names(known_deltas), drop = FALSE] %*% known_deltas
  I_minus_B_UU <- diag(length(unknown_genes)) - B_UU
  imputed <- base::solve(I_minus_B_UU, driving_force)
  return(imputed[,1])
}


#' Impute Absolute Expressions
#' Solves the linear system in the absolute space: x = α + Bx.
#' @noRd
.impute_absolute <- function(B_full, known_expressions) {
  # B_full is targets x (sources + Intercept).
  genes <- rownames(B_full)

  unknown_genes <- setdiff(genes, names(known_expressions))
  if (length(unknown_genes) == 0) return(numeric(0))

  # B_no_intercept should be the square, gene-only part of the matrix.
  B_no_intercept <- B_full[genes, genes, drop = FALSE]
  B_UU <- B_no_intercept[unknown_genes, unknown_genes, drop = FALSE]

  epsilon <- 1e-8
  eigen_vals <- base::eigen(B_UU, only.values = TRUE)$values
  if (any(abs(eigen_vals) >= 1 - epsilon)) {
    warning("System may be unstable (max |eigenvalue| >= 1 - epsilon).")
  }

  known_with_intercept <- c(1, known_expressions)
  names(known_with_intercept) <- c('Intercept', names(known_expressions))

  # The driving force includes the intercept term.
  driving_force <- B_full[unknown_genes, names(known_with_intercept), drop = FALSE] %*% known_with_intercept
  I_minus_B_UU <- diag(length(unknown_genes)) - B_UU
  imputed <- base::solve(I_minus_B_UU, driving_force)
  return(imputed[,1])
}


#' Predict Standard Perturbation Effect
#'
#' Predicts the perturbation effect for a set of KOs by propagating the
#' initial perturbation through the network structure. Allows for modeling in
#' either the "delta" or "absolute" expression space.
#'
#' @param B A numeric matrix of regulatory coefficients (sources x targets).
#'   This is the output of \code{\link{fit_expression_model}}.
#' @param ko_expressions A named vector of expression levels for perturbed genes.
#' @param wt_expressions A named vector of wild-type expression levels.
#' @param max_dist Maximum distance for effect propagation. Default is `Inf`.
#'   Distances are calculated on the functional graph derived from B.
#' @param space The modeling space. delta: propagate in expression change space. absolute: propagate in expression space.
#' @examples
#' # --- 0. SETUP: Load Libraries & Define Ground Truth ---
#' library(dplyr)
#' library(igraph)
#' library(CausalGRN)
#'
#' # Define all simulation parameters upfront
#' a <- b <- 1
#' sd <- 2
#' s <- -4
#' nwt <- 1e5
#' npt <- 1e4
#' ko_efficacy <- 0.9 # 90% knockdown efficiency
#'
#' # --- 1. DATA SIMULATION: Create the Toy Example ---
#' set.seed(123)
#'
#' # Generate Wild-Type (WT) Data
#' x_latent_wt <- rnorm(nwt, 0, sd)
#' y_latent_wt <- rnorm(nwt, a * x_latent_wt + s, sd)
#' z_latent_wt <- rnorm(nwt, b * (y_latent_wt - s), sd)
#' wt_counts <- cbind(
#'   A = rpois(nwt, exp(x_latent_wt)),
#'   B = rpois(nwt, exp(y_latent_wt)),
#'   C = rpois(nwt, exp(z_latent_wt))
#' )
#'
#' # Generate Perturb-A Data
#' x_latent_koA <- rnorm(npt, 0, sd) + log(1 - ko_efficacy)
#' y_latent_koA <- rnorm(npt, a * x_latent_koA + s, sd)
#' z_latent_koA <- rnorm(npt, b * (y_latent_koA - s), sd)
#' koA_counts <- cbind(
#'   A = rpois(npt, exp(x_latent_koA)),
#'   B = rpois(npt, exp(y_latent_koA)),
#'   C = rpois(npt, exp(z_latent_koA))
#' )
#'
#' # Generate Perturb-B Data
#' x_latent_koB <- rnorm(npt, 0, sd)
#' y_latent_koB <- rnorm(npt, a * x_latent_koB + s, sd) + log(1 - ko_efficacy)
#' z_latent_koB <- rnorm(npt, b * (y_latent_koB - s), sd)
#' koB_counts <- cbind(
#'   A = rpois(npt, exp(x_latent_koB)),
#'   B = rpois(npt, exp(y_latent_koB)),
#'   C = rpois(npt, exp(z_latent_koB))
#' )
#'
#' # Generate Perturb-C Data
#' x_latent_koC <- rnorm(npt, 0, sd)
#' y_latent_koC <- rnorm(npt, a * x_latent_koC + s, sd)
#' z_latent_koC <- rnorm(npt, b * (y_latent_koC - s), sd) + log(1 - ko_efficacy)
#' koC_counts <- cbind(
#'   A = rpois(npt, exp(x_latent_koC)),
#'   B = rpois(npt, exp(y_latent_koC)),
#'   C = rpois(npt, exp(z_latent_koC))
#' )
#'
#' # --- 2. PREPARE INPUTS ---
#' count <- rbind(wt_counts, koA_counts, koB_counts, koC_counts)
#' group <- c(rep('WT', nwt), rep('A', npt), rep('B', npt), rep('C', npt))
#' Y <- scale(log1p(count), center = TRUE, scale = TRUE)
#' colnames(count) <- colnames(Y) <- c('A', 'B', 'C')
#'
#' # --- 3. Fit model (using WT and KO A data) ---
#' train_idx <- which(group %in% c('WT', 'A'))
#' skel <- infer_skeleton(
#'   count[train_idx, ],
#'   Y[train_idx, ],
#'   alpha = 0.05,
#'   min_abspcor = 0,
#'   ncores = 1
#' )
#' stat <- calc_perturbation_effect(Y[train_idx, ], group[train_idx], ncores = 1)
#' causal_graph <- infer_causalgrn(skel$graph, stat, alpha = 0.05)
#' B <- fit_expression_model(
#'   Y[train_idx, ],
#'   group[train_idx],
#'   graph = causal_graph,
#'   ncores = 1,
#'   method = 'lm'
#' )
#'
#' # --- 4. Predict effects for B and C knockout ---
#' wt_expressions <- colMeans(Y[group == 'WT', ])
#' ko_expressions <- c(
#'   'B' = mean(Y[group == 'B', 'B']),
#'   'C' = mean(Y[group == 'C', 'C'])
#' )
#' pred_effects <- predict_standard_effect(B, ko_expressions, wt_expressions)
#' print(pred_effects)
#' @return A numeric matrix of predicted delta values (KOs x genes).
#' @export
predict_standard_effect <- function(B, ko_expressions, wt_expressions, max_dist = Inf, space = c("delta", "absolute")) {
  # --- 1. Input Validation and Setup ---
  space <- match.arg(space)
  genes <- names(wt_expressions)
  ko_genes <- names(ko_expressions)
  stopifnot(
    is.matrix(B), is.numeric(wt_expressions), is.numeric(ko_expressions),
    identical(colnames(B), genes),
    identical(rownames(B), c('Intercept', genes)),
    all(ko_genes %in% genes)
  )

  # --- 2. Prepare B Matrix and Functional Graph ---
  B_propagator <- t(B[genes, , drop = FALSE])
  adj_matrix_functional <- t(B_propagator != 0)
  graph_for_distances <- igraph::graph_from_adjacency_matrix(adj_matrix_functional, mode = 'directed')

  # --- 3. Main Prediction Loop ---
  pred_delta_list <- lapply(ko_genes, function(ko_gene) {

    distances <- igraph::distances(graph_for_distances, v = ko_gene, mode = 'out')[1, ]
    unchanged_genes <- names(which(is.infinite(distances) | distances > max_dist))

    if (space == "delta") {
      known_deltas <- c()
      known_deltas[ko_gene] <- ko_expressions[ko_gene] - wt_expressions[ko_gene]
      known_deltas[unchanged_genes] <- 0

      imputed_deltas <- .impute_deltas(B_propagator = B_propagator, known_deltas = known_deltas)
      if (length(imputed_deltas) > 0) known_deltas[names(imputed_deltas)] <- imputed_deltas
      stopifnot(
        length(known_deltas) == length(genes),
        setequal(names(known_deltas), genes)
      )
      return(known_deltas[genes])

    } else if (space == 'absolute') {
      pred_ko_abs <- c()
      pred_ko_abs[ko_gene] <- ko_expressions[ko_gene]
      pred_ko_abs[unchanged_genes] <- wt_expressions[unchanged_genes]

      imputed_abs <- .impute_absolute(B_full = t(B), known_expressions = pred_ko_abs)
      if (length(imputed_abs) > 0) pred_ko_abs[names(imputed_abs)] <- imputed_abs

      # Convert the final absolute prediction to a delta
      stopifnot(
        length(pred_ko_abs) == length(genes),
        setequal(names(pred_ko_abs), genes)
      )
      return(pred_ko_abs[genes] - wt_expressions)
    }
  })

  # --- 4. Assemble and Return Final Matrix ---
  pred_matrix <- do.call(rbind, pred_delta_list)
  rownames(pred_matrix) <- ko_genes

  return(pred_matrix)
}
