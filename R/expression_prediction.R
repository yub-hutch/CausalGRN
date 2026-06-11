#' Fit Gene Expression Model
#'
#' Fits a predictive model for each gene's expression based on a given graph
#' structure and regression method.
#'
#' @param Y A numeric matrix of expression data (cells x genes).
#' @param group Named character vector of cell labels: 'WT' for wild-type cells,
#'   perturbed gene for perturbed cell.
#' @param graph An igraph object or the character string 'all' specifying the
#'   regulatory structure.
#' @param ncores The number of cores for parallel computation.
#' @param method The regression method to use. One of 'lm', 'lasso', or 'ridge'.
#' @return A numeric matrix `B` where rows are sources ('Intercept' and genes)
#'   and columns are targets (genes).
#' @export
fit_expression_model <- function(
    Y, group, graph, ncores, method = c('lm', 'lasso', 'ridge')
) {
  .check_perturbation_effect_inputs(Y = Y, group = group)
  .check_ncores(ncores)

  p <- ncol(Y)
  genes <- colnames(Y)
  method <- match.arg(method)
  if (method != 'lm' && !requireNamespace("glmnet", quietly = TRUE)) {
    stop(
      "Package 'glmnet' must be installed when method is 'lasso' or 'ridge'.",
      call. = FALSE
    )
  }

  if (inherits(graph, "igraph")) {
    .check_igraph(graph, nodes = genes)
  } else if (!identical(graph, 'all')) {
    stop("'graph' must be an igraph object or the string 'all'.", call. = FALSE)
  }

  message(
    'Fitting models for ', p, ' genes with ', nrow(Y),
    ' cells using method: ', method, ' ...'
  )

  if (identical(graph, 'all')) {
    adj_matrix <- matrix(TRUE, p, p, dimnames = list(genes, genes))
    diag(adj_matrix) <- FALSE
    graph <- igraph::graph_from_adjacency_matrix(adj_matrix, mode = 'directed')
  }

  model_list <- .parallel_lapply(
    genes,
    function(gene) {
      coef_vector <- setNames(rep(0, p + 1), c('Intercept', genes))

      predictors <- igraph::neighbors(graph, gene, mode = 'in')$name
      samples <- which(group != gene)

      if (length(predictors) == 0L) {
        coef_vector['Intercept'] <- mean(Y[samples, gene])
        return(coef_vector)
      }

      if (length(predictors) == 1L || method == 'lm') {
        fit <- stats::lm(Y[samples, gene] ~ Y[samples, predictors])
        fit_coefs <- stats::coef(fit)

        if (anyNA(fit_coefs)) {
          warning(
            "lm failed for gene ", gene, " - fitting intercept only.",
            call. = FALSE
          )
          coef_vector['Intercept'] <- mean(Y[samples, gene])
          return(coef_vector)
        }

        coef_vector[c('Intercept', predictors)] <- fit_coefs
        return(coef_vector)
      }

      alpha_val <- if (method == 'lasso') 1 else 0
      cvfit <- glmnet::cv.glmnet(
        x = Y[samples, predictors, drop = FALSE],
        y = Y[samples, gene],
        family = 'gaussian',
        nfolds = 5,
        alpha = alpha_val,
        intercept = TRUE
      )
      fit_coefs <- as.matrix(stats::coef(cvfit, s = 'lambda.min'))

      coef_vector['Intercept'] <- fit_coefs[1, 1]
      coef_vector[predictors] <- fit_coefs[-1, 1]
      return(coef_vector)
    },
    ncores = ncores,
    export = c('Y', 'group', 'graph', 'p', 'genes', 'method')
  )

  B <- do.call(cbind, model_list)
  colnames(B) <- genes

  if (any(is.na(B))) {
    stop("Fitted coefficient matrix contains missing values.", call. = FALSE)
  }
  return(B)
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
  if (length(unknown_genes) == 0) {
    return(numeric(0))
  }

  B_UU <- B_propagator[unknown_genes, unknown_genes, drop = FALSE]

  epsilon <- 1e-8
  eigen_vals <- base::eigen(B_UU, only.values = TRUE)$values
  if (any(abs(eigen_vals) >= 1 - epsilon)) {
    warning("System may be unstable (max |eigenvalue| >= 1 - epsilon).")
  }

  driving_force <- B_propagator[
    unknown_genes,
    names(known_deltas),
    drop = FALSE
  ] %*% known_deltas
  I_minus_B_UU <- diag(length(unknown_genes)) - B_UU
  imputed <- base::solve(I_minus_B_UU, driving_force)
  return(imputed[, 1])
}


#' Predict Standard Perturbation Effect
#'
#' Predicts the perturbation effect for a set of KOs by propagating the
#' initial perturbation through the network structure in delta space.
#'
#' @param B A numeric matrix of regulatory coefficients (sources x targets).
#'   This is the output of \code{\link{fit_expression_model}}.
#' @param ko_expressions A named vector of expression levels for perturbed genes.
#' @param wt_expressions A named vector of wild-type expression levels.
#' @param max_dist Maximum distance for effect propagation. Default is `Inf`.
#'   Distances are calculated on the functional graph derived from B.
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
#' rownames(count) <- paste0("cell", seq_len(nrow(count)))
#' names(group) <- rownames(count)
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
#' pred_effects <- predict_perturbation_effect(B, ko_expressions, wt_expressions)
#' print(pred_effects)
#' @return A numeric matrix of predicted delta values (KOs x genes).
#' @export
predict_perturbation_effect <- function(
    B, ko_expressions, wt_expressions, max_dist = Inf
) {
  .check_expression_vector(wt_expressions, arg = 'wt_expressions')
  .check_expression_vector(ko_expressions, arg = 'ko_expressions')

  genes <- names(wt_expressions)
  ko_genes <- names(ko_expressions)

  .check_expression_model_matrix(B, genes = genes)
  if (!all(ko_genes %in% genes)) {
    stop(
      "'ko_expressions' names must be included in 'wt_expressions'.",
      call. = FALSE
    )
  }
  if (
    length(max_dist) != 1L || !is.numeric(max_dist) ||
      is.na(max_dist) || max_dist < 0
  ) {
    stop("'max_dist' must be a single non-negative number or Inf.", call. = FALSE)
  }

  B_propagator <- t(B[genes, , drop = FALSE])
  adj_matrix_functional <- t(B_propagator != 0)
  graph_for_distances <- igraph::graph_from_adjacency_matrix(
    adj_matrix_functional,
    mode = 'directed'
  )

  pred_delta_list <- lapply(ko_genes, function(ko_gene) {
    distances <- igraph::distances(
      graph_for_distances,
      v = ko_gene,
      mode = 'out'
    )[1, ]
    unchanged_genes <- names(which(
      is.infinite(distances) | distances > max_dist
    ))

    known_deltas <- numeric(0)
    known_deltas[ko_gene] <- ko_expressions[ko_gene] - wt_expressions[ko_gene]
    known_deltas[unchanged_genes] <- 0

    imputed_deltas <- .impute_deltas(
      B_propagator = B_propagator,
      known_deltas = known_deltas
    )
    if (length(imputed_deltas) > 0) {
      known_deltas[names(imputed_deltas)] <- imputed_deltas
    }
    if (
      length(known_deltas) != length(genes) ||
        !setequal(names(known_deltas), genes)
    ) {
      stop("Prediction did not return a value for all genes.", call. = FALSE)
    }
    return(known_deltas[genes])
  })

  pred_matrix <- do.call(rbind, pred_delta_list)
  rownames(pred_matrix) <- ko_genes

  return(pred_matrix)
}
