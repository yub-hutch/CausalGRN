#' Fit a Two-Stage Portability Prediction Model
#'
#' Trains two model components: a standard GRN model (B_matrix) for the target
#' dataset, and a flexible set of portability models to predict the SCOPE
#' effect using information from source datasets.
#'
#' @param Y A numeric matrix of expression data for the target dataset.
#' @param group A character vector for the cell groups in Y.
#' @param graph An igraph object for the target dataset.
#' @param source_effects A named list of data frames with KO effects from source datasets.
#' @param pname The name of the hub node (e.g., "PC1").
#' @param pgenes A character vector of the signature genes.
#' @param ncores The number of cores for parallel computation.
#' @return A list containing the `B_matrix` and a list of trained `portability_models`.
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

  message("--- Fitting Two-Stage Portability Model ---")

  # --- 2. Learn the B Matrix for the Target Dataset ---
  message("Step 1/2: Learning main B matrix for the target dataset...")

  B_matrix <- fit_expression_model(
    Y = Y, group = group, graph = graph,
    ncores = ncores,
    method = 'lm'
  )

  # --- 3. Learn the Portability Models (Stage 1 Regressors) ---
  message("Step 2/2: Learning flexible portability models (Source -> Target)...")

  wt_mean_pname <- mean(Y[group == "WT", pname])
  target_effects <- dplyr::tibble(
    ko = setdiff(unique(group), "WT"),
    target_delta = sapply(ko, \(k) mean(Y[group == k, pname]) - wt_mean_pname)
  )

  source_names <- names(source_effects)
  portability_models <- list()

  # a. Fit a model using all available source predictors
  all_dfs_to_join <- c(list(target_effects), source_effects)
  training_df_both <- Reduce(\(x, y) dplyr::inner_join(x, y, by = "ko"), all_dfs_to_join)

  message(paste("  - Fitting 'both' model using", nrow(training_df_both), "complete KO observations."))

  stopifnot(
    "Not enough complete KO observations to fit 'both' model. Need more KOs than predictors." =
      nrow(training_df_both) > length(source_names)
  )
  formula_both <- as.formula(paste("target_delta ~", paste(source_names, collapse = " + ")))
  portability_models[['both']] <- stats::lm(formula_both, data = training_df_both)

  # b. Fit a model for each single source predictor
  for (s_name in source_names) {
    training_df_single <- dplyr::inner_join(target_effects, source_effects[[s_name]], by = "ko")

    message(paste("  - Fitting '", s_name, "' model using", nrow(training_df_single), "complete KO observations."))

    stopifnot(
      "Not enough complete KO observations to fit single-source model. Need at least 2." =
        nrow(training_df_single) > 1
    )
    formula_single <- as.formula(paste("target_delta ~", s_name))
    portability_models[[s_name]] <- stats::lm(formula_single, data = training_df_single)
  }

  message("--- Component training complete. ---")

  # --- 4. Return Both Components ---
  return(list(
    B_matrix = B_matrix,
    portability_models = portability_models
  ))
}


#' Predict Perturbation Effect Using a Two-Stage Portability Model
#'
#' @param B A numeric matrix of regulatory coefficients.
#' @param portability_models A named list of trained `lm` objects.
#' @param source_effects A named list of data frames with predictor variables.
#' @param ko_expressions A named vector of expression levels for the test KOs.
#' @param wt_expressions A named vector of WT expressions for the target dataset.
#' @param pname The name of the hub node.
#' @param max_dist Maximum distance for effect propagation.
#' @param ncores Number of cores for parallel computation.
#' @return A numeric matrix of predicted delta values for the test KOs.
#' @export
predict_standard_effect_with_gp <- function(
    B, portability_models, source_effects,
    ko_expressions, wt_expressions, pname,
    max_dist = Inf, ncores = 1
) {
  # --- 1. Setup ---
  genes <- names(wt_expressions)
  ko_genes <- names(ko_expressions)
  stopifnot(all(ko_genes %in% source_effects[[1]]$ko)) # A basic check

  source_names <- names(portability_models)[names(portability_models) != 'both']

  # --- 2. Stage 1: Predict the SCOPE Effect for Test KOs ---
  message("  Stage 1: Predicting SCOPE effect for each test KO...")
  predicted_pname_deltas <- sapply(ko_genes, function(ko) {

    available_sources <- c()
    for (s_name in source_names) {
      if (ko %in% source_effects[[s_name]]$ko && !is.na(source_effects[[s_name]][source_effects[[s_name]]$ko == ko, s_name])) {
        available_sources <- c(available_sources, s_name)
      }
    }

    if (length(available_sources) == 0) return(0)

    ko_data_list <- lapply(available_sources, function(s_name) {
      source_effects[[s_name]] %>% dplyr::filter(ko == !!ko)
    })
    ko_data <- Reduce(\(x, y) dplyr::full_join(x, y, by = "ko"), ko_data_list)

    if (length(available_sources) == length(source_names) && 'both' %in% names(portability_models)) {
      return(stats::predict(portability_models[['both']], newdata = ko_data))
    } else {
      first_available_model <- intersect(available_sources, names(portability_models))[1]
      if (!is.na(first_available_model)) {
        return(stats::predict(portability_models[[first_available_model]], newdata = ko_data))
      }
    }
    return(0)
  })
  names(predicted_pname_deltas) <- ko_genes

  # --- 3. Stage 2: Propagate Effects through the GRN ---
  message("  Stage 2: Propagating effects through the GRN...")
  B_propagator <- t(B[genes, , drop = FALSE])

  pred_delta_list <- pbmcapply::pbmclapply(ko_genes, function(ko_gene) {

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

    # Return in the correct, consistent order
    return(known_deltas[genes])

  }, mc.cores = ncores)

  # --- 4. Assemble and Return Final Matrix ---
  pred_matrix <- do.call(rbind, pred_delta_list)
  rownames(pred_matrix) <- ko_genes

  return(pred_matrix)
}
