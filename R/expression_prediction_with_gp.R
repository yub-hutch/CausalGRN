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
#' @param constrain_graph Logical. If TRUE, constrain pgenes to only have pname as a parent before fitting the GRN.
#' @param alpha Significance level threshold (F-statistic p-value) to determine if a portability model is considered valid.
#' @return A list containing `B_matrix` and `portability_models`. `portability_models` contains a model for each combination of sources, plus fallback values (`fallback_pname_delta`, `fallback_pgenes_deltas`).
#' @export
fit_expression_model_with_gp <- function(
    Y, group, graph, source_effects, pname, pgenes, ncores,
    constrain_graph = TRUE,
    alpha = 0.05
) {

  all_nodes <- colnames(Y)
  stopifnot(
    is.list(source_effects) && !is.data.frame(source_effects),
    inherits(graph, "igraph"),
    setequal(all_nodes, igraph::V(graph)$name),
    !igraph::any_loop(graph),
    all(c(pname, pgenes) %in% all_nodes)
  )

  # Fit the B Matrix for the Target Dataset
  fit_graph <- if (constrain_graph) {
    adj_matrix <- as.matrix(igraph::as_adjacency_matrix(graph))
    adj_matrix[setdiff(rownames(adj_matrix), pname), pgenes] <- 0
    igraph::graph_from_adjacency_matrix(adj_matrix, mode = 'directed')
  } else {
    graph
  }

  B_matrix <- fit_expression_model(
    Y = Y, group = group, graph = fit_graph,
    method = 'lm', ncores = ncores
  )

  # Learn the Portability Models (Stage 1 Regressors)
  training_kos <- setdiff(unique(group), "WT")
  wt_means_all <- colMeans(Y[group == "WT", , drop = FALSE])

  target_effects <- dplyr::tibble(
    ko = training_kos,
    target_delta = sapply(training_kos, function(k) {
      mean(Y[group == k, pname]) - wt_means_all[pname]
    })
  )

  source_names <- names(source_effects)
  portability_models <- list()

  # Calculate and store fallback values
  portability_models[['fallback_pname_delta']] <- mean(
    target_effects$target_delta, na.rm = TRUE
  )

  # Calculate mean deltas for all genes, but only store those for pgenes
  ko_deltas_all_list <- lapply(training_kos, function(k) {
    colMeans(Y[group == k, , drop = FALSE]) - wt_means_all
  })
  ko_deltas_all_matrix <- do.call(rbind, ko_deltas_all_list)
  mean_all_gene_deltas <- colMeans(ko_deltas_all_matrix, na.rm = TRUE)

  portability_models[['fallback_pgenes_deltas']] <- mean_all_gene_deltas[pgenes]

  # Create master df starting from target KOs and left-joining source data
  all_dfs_to_join <- c(list(target_effects), source_effects)
  training_df_master <- Reduce(
    function(x, y) dplyr::left_join(x, y, by = "ko"),
    all_dfs_to_join
  )

  # Generate all 2^N - 1 combinations
  all_combos_list <- lapply(1:length(source_names), function(k) {
    utils::combn(source_names, k, simplify = FALSE)
  })
  all_combos <- unlist(all_combos_list, recursive = FALSE)

  for (combo in all_combos) {
    present_sources <- combo
    model_name <- paste(sort(present_sources), collapse = "_and_")

    df_combo <- tidyr::drop_na(training_df_master, dplyr::all_of(present_sources))

    if (nrow(df_combo) > length(present_sources)) {
      formula_combo <- stats::as.formula(
        paste("target_delta ~", paste(present_sources, collapse = " + "))
      )
      lm_combo <- stats::lm(formula_combo, data = df_combo)

      s <- summary(lm_combo)
      fstat <- s$fstatistic
      p_value <- stats::pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE)
      is_significant <- !is.na(p_value) && p_value < alpha

      obs_scale <- stats::sd(df_combo$target_delta)
      pred_scale <- stats::sd(stats::predict(lm_combo, newdata = df_combo))
      scale_factor <- ifelse(pred_scale > 1e-6, obs_scale / pred_scale, 1.0)

      portability_models[[model_name]] <- list(
        model = lm_combo,
        scale_factor = scale_factor,
        is_significant = is_significant
      )
    } else {
      portability_models[[model_name]] <- list(
        model = NULL,
        scale_factor = 1.0,
        is_significant = FALSE
      )
    }
  }

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
#' @param portability_models A named list of trained models and fallback values.
#' @param source_effects A named list of data frames with predictor variables.
#' @param ko_expressions A named vector of expression levels for the test KOs.
#' @param wt_expressions A named vector of WT expressions for the target dataset.
#' @param pname The name of the hub node.
#' @param pgenes A character vector of the signature genes.
#' @param graph An igraph object, required if `ancestors_only = TRUE`.
#' @param scale_pname Logical. If TRUE, apply the stored variance calibration.
#' @param ancestors_only Logical. If TRUE, only apply portability models to KOs that are ancestors of `pname` in the `graph`.
#' @param fallback_strategy Character. Action if no source data is found or KO is not an ancestor: "pname_mean_propagate", "pname_zero_propagate", or "pgenes_mean_propagate".
#' @return A numeric matrix of predicted delta values for the test KOs.
#' @export
predict_standard_effect_with_gp <- function(
    B, portability_models, source_effects,
    ko_expressions, wt_expressions, pname, pgenes,
    graph = NULL,
    scale_pname = TRUE,
    ancestors_only = FALSE,
    fallback_strategy = "pname_mean_propagate"
) {

  genes <- names(wt_expressions)
  ko_genes <- names(ko_expressions)
  all_source_names <- names(source_effects)
  B_propagator <- t(B[genes, , drop = FALSE])

  # --- 1. Validate Arguments and Get Fallback Values ---
  fallback_strategy <- match.arg(
    fallback_strategy,
    c("pname_mean_propagate", "pname_zero_propagate", "pgenes_mean_propagate")
  )

  fallback_pname_delta <- portability_models[['fallback_pname_delta']]
  if (is.null(fallback_pname_delta)) {
    stop("Corrupted 'portability_models': 'fallback_pname_delta' is missing.")
  }

  fallback_pgenes_deltas <- portability_models[['fallback_pgenes_deltas']]
  if (fallback_strategy == "pgenes_mean_propagate") {
    if (is.null(fallback_pgenes_deltas)) {
      stop("Corrupted 'portability_models': 'fallback_pgenes_deltas' is missing.")
    }
    if (!setequal(pgenes, names(fallback_pgenes_deltas))) {
      stop("'pgenes' do not match the keys in 'fallback_pgenes_deltas'.")
    }
  }

  # --- 2. Pre-calculate Ancestor Status (More Efficient) ---
  ancestor_map <- NULL
  if (ancestors_only) {
    if (!inherits(graph, "igraph")) {
      stop("'graph' must be a valid igraph object when 'ancestors_only = TRUE'.")
    }
    dists <- igraph::distances(graph, v = ko_genes, to = pname, mode = "out")
    ancestor_map <- setNames(is.finite(dists[, 1]), ko_genes)
  }

  # --- 3. Loop Through KOs, Predict, and Propagate ---
  pred_delta_list <- lapply(ko_genes, function(ko_gene) {

    use_fallback <- FALSE
    predicted_pname_delta <- NA_real_

    # Check ancestor condition if specified
    if (ancestors_only && !ancestor_map[[ko_gene]]) {
      use_fallback <- TRUE
    }

    # If not falling back yet, try to find a portability model
    if (!use_fallback) {
      available_sources <- c()
      for (s_name in all_source_names) {
        if (ko_gene %in% source_effects[[s_name]]$ko && !is.na(source_effects[[s_name]][source_effects[[s_name]]$ko == ko_gene, s_name])) {
          available_sources <- c(available_sources, s_name)
        }
      }

      if (length(available_sources) == 0) {
        use_fallback <- TRUE
      } else {
        model_name <- paste(sort(available_sources), collapse = "_and_")
        model_entry <- portability_models[[model_name]]

        if (is.null(model_entry) || !isTRUE(model_entry$is_significant)) {
          use_fallback <- TRUE
        } else {
          ko_data_list <- lapply(available_sources, function(s_name) {
            source_effects[[s_name]] %>% dplyr::filter(ko == !!ko_gene)
          })
          ko_data <- Reduce(function(x, y) dplyr::full_join(x, y, by = "ko"), ko_data_list)

          raw_pred <- stats::predict(model_entry$model, newdata = ko_data)

          predicted_pname_delta <- if (scale_pname) {
            raw_pred * model_entry$scale_factor
          } else {
            raw_pred
          }
        }
      }
    }

    # Apply fallback logic to set the initial predicted_pname_delta
    if (use_fallback) {
      if (fallback_strategy == "pname_mean_propagate" || fallback_strategy == "pgenes_mean_propagate") {
        predicted_pname_delta <- fallback_pname_delta
      } else { # "pname_zero_propagate"
        predicted_pname_delta <- 0.0
      }
    }

    # Propagate effects through the GRN
    known_deltas <- c()
    known_deltas[ko_gene] <- ko_expressions[ko_gene] - wt_expressions[ko_gene]
    known_deltas[pname] <- predicted_pname_delta

    # If using pgenes_mean_propagate, set those values *before* imputation
    if (use_fallback && fallback_strategy == "pgenes_mean_propagate") {
      # Set pgenes, but NEVER overwrite the ko_gene's own delta
      pgenes_to_set <- setdiff(pgenes, ko_gene)
      if (length(pgenes_to_set) > 0) {
        known_deltas[pgenes_to_set] <- fallback_pgenes_deltas[pgenes_to_set]
      }
    }

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
