#' Fit a Two-Stage Portability Prediction Model
#'
#' Trains a GRN model and a set of hierarchically inclusive portability models
#' for all combinations of source datasets.
#'
#' @param Y A numeric matrix of expression data for the target dataset.
#' @param group A character vector for the cell groups in Y.
#' @param graph An igraph object for the target dataset.
#' @param source_effects A named list of data frames with KO effects from source datasets.
#' @param pname The name of the hub node (e.g., "PC1").
#' @param pgenes A character vector of the signature genes.
#' @param ncores The number of cores for parallel computation.
#' @return A list containing `B_matrix` and `portability_models`.
#'   `portability_models` contains a model for each combination of sources,
#'   plus a `mean_fallback` value.
#' @export
fit_expression_model_with_gp <- function(
    Y, group, graph, source_effects, pname, pgenes, ncores
) {

  # --- 1. Parameter Validation ---
  all_nodes <- colnames(Y)
  stopifnot(is.list(source_effects) && !is.data.frame(source_effects))
  stopifnot(
    inherits(graph, "igraph"),
    setequal(all_nodes, igraph::V(graph)$name),
    !igraph::any_loop(graph),
    all(c(pname, pgenes) %in% all_nodes)
  )

  # --- 2. Learn the B Matrix for the Target Dataset ---
  adj_matrix <- as.matrix(igraph::as_adjacency_matrix(graph))
  adj_matrix[setdiff(rownames(adj_matrix), pname), pgenes] <- 0
  constrained_graph <- igraph::graph_from_adjacency_matrix(adj_matrix, mode = 'directed')

  B_matrix <- fit_expression_model(
    Y = Y, group = group, graph = constrained_graph,
    method = 'lm', ncores = ncores
  )

  # --- 3. Learn the Portability Models (Stage 1 Regressors) ---
  wt_mean_pname <- mean(Y[group == "WT", pname])
  target_effects <- dplyr::tibble(
    ko = setdiff(unique(group), "WT"),
    target_delta = sapply(setdiff(unique(group), "WT"), function(k) {
      mean(Y[group == k, pname]) - wt_mean_pname
    })
  )

  source_names <- names(source_effects)
  portability_models <- list()

  # Calculate and store the mean_fallback for KOs with no source data
  mean_fallback <- mean(target_effects$target_delta, na.rm = TRUE)
  portability_models[['mean_fallback']] <- mean_fallback

  # Create master df starting from target KOs and left-joining source data
  all_dfs_to_join <- c(list(target_effects), source_effects)
  training_df_master <- Reduce(
    function(x, y) dplyr::left_join(x, y, by = "ko"),
    all_dfs_to_join
  )

  # Generate all 2^N - 1 combinations (power set, excluding empty set)
  all_combos_list <- lapply(1:length(source_names), function(k) {
    utils::combn(source_names, k, simplify = FALSE)
  })
  all_combos <- unlist(all_combos_list, recursive = FALSE)

  # Loop through every combination and fit a model
  for (combo in all_combos) {

    present_sources <- combo
    model_name <- paste(sort(present_sources), collapse = "_and_")

    # Filter for KOs that have *at least* the present_sources
    df_combo <- tidyr::drop_na(training_df_master, dplyr::all_of(present_sources))

    # Check if we have enough data to fit this model (more KOs than predictors)
    if (nrow(df_combo) > length(present_sources)) {

      formula_combo <- stats::as.formula(paste("target_delta ~", paste(present_sources, collapse = " + ")))

      lm_combo <- stats::lm(formula_combo, data = df_combo)

      # Calculate variance calibration factor
      obs_scale <- stats::sd(df_combo$target_delta)
      pred_scale <- stats::sd(stats::predict(lm_combo, newdata = df_combo))
      scale_factor <- ifelse(pred_scale > 1e-6, obs_scale / pred_scale, 1.0)

      # Store the calibrated model
      portability_models[[model_name]] <- list(
        model = lm_combo,
        scale_factor = scale_factor
      )

    } else {
      stop(paste(
        "Not enough observations to fit portability model '", model_name, "'.",
        "Need >", length(present_sources), "KOs, but found", nrow(df_combo)
      ))
    }
  }

  # --- 4. Return Both Components ---
  return(list(
    B_matrix = B_matrix,
    portability_models = portability_models
  ))
}


#' Predict Perturbation Effect Using a Two-Stage Portability Model
#'
#' Predicts the perturbation effect by matching test KOs to the best-available
#' pre-trained portability model and propagating the effect.
#'
#' @param B A numeric matrix of regulatory coefficients.
#' @param portability_models A named list of trained `lm` objects and a fallback.
#' @param source_effects A named list of data frames with predictor variables.
#' @param ko_expressions A named vector of expression levels for the test KOs.
#' @param wt_expressions A named vector of WT expressions for the target dataset.
#' @param pname The name of the hub node.
#' @return A numeric matrix of predicted delta values for the test KOs.
#' @export
predict_standard_effect_with_gp <- function(
    B, portability_models, source_effects,
    ko_expressions, wt_expressions, pname
) {
  # --- 1. Setup ---
  genes <- names(wt_expressions)
  ko_genes <- names(ko_expressions)

  all_source_names <- names(source_effects)

  # Get the mean_fallback value, stopping if it's missing.
  mean_fallback <- portability_models[['mean_fallback']]
  if (is.null(mean_fallback)) {
    stop("Corrupted 'portability_models': 'mean_fallback' value is missing.")
  }

  # --- 2. Stage 1: Predict the SCOPE Effect for Test KOs ---
  predicted_pname_deltas <- sapply(ko_genes, function(ko) {

    # Find all available sources for this specific KO
    available_sources <- c()
    for (s_name in all_source_names) {
      if (ko %in% source_effects[[s_name]]$ko && !is.na(source_effects[[s_name]][source_effects[[s_name]]$ko == ko, s_name])) {
        available_sources <- c(available_sources, s_name)
      }
    }

    # If no sources are available, use the training mean
    if (length(available_sources) == 0) {
      return(mean_fallback)
    }

    # Build the exact model name based on the available sources
    model_name <- paste(sort(available_sources), collapse = "_and_")

    # The 'fit' function should have created this model. If not, it's an error.
    stopifnot(
      "Internal error: No matching portability model found." =
        model_name %in% names(portability_models)
    )

    # The exact model exists, so we use it
    ko_data_list <- lapply(available_sources, function(s_name) {
      source_effects[[s_name]] %>% dplyr::filter(ko == !!ko)
    })
    ko_data <- Reduce(
      function(x, y) dplyr::full_join(x, y, by = "ko"),
      ko_data_list
    )

    # Get the calibrated prediction
    model_fit <- portability_models[[model_name]]$model
    scale_factor <- portability_models[[model_name]]$scale_factor
    raw_pred <- stats::predict(model_fit, newdata = ko_data)

    return(raw_pred * scale_factor)
  })
  names(predicted_pname_deltas) <- ko_genes

  # --- 3. Stage 2: Propagate Effects through the GRN ---
  B_propagator <- t(B[genes, , drop = FALSE])

  pred_delta_list <- lapply(ko_genes, function(ko_gene) {

    known_deltas <- c()
    known_deltas[ko_gene] <- ko_expressions[ko_gene] - wt_expressions[ko_gene]
    known_deltas[pname] <- predicted_pname_deltas[ko_gene]

    imputed_deltas <- .impute_deltas(B_propagator = B_propagator, known_deltas = known_deltas)
    if (length(imputed_deltas) > 0) known_deltas[names(imputed_deltas)] <- imputed_deltas

    stopifnot(
      "Prediction did not return a value for all genes." =
        setequal(names(known_deltas), genes),
      "Prediction resulted in NA values." =
        !any(is.na(known_deltas))
    )

    return(known_deltas[genes])
  })

  # --- 4. Assemble and Return Final Matrix ---
  pred_matrix <- do.call(rbind, pred_delta_list)
  rownames(pred_matrix) <- ko_genes

  return(pred_matrix)
}
