#' Fit Expression Model With Gene Program
#'
#' Trains a GRN model and a set of hierarchically inclusive portability models
#' for all combinations of source datasets.
#'
#' @param Y A numeric matrix of expression data for the target dataset.
#' @param group Named character vector for the cell groups in Y.
#' @param graph An igraph object for the target dataset.
#' @param source_effects A named list of source-effect data frames. Each list
#'   name is a source name. Each data frame must contain exactly two columns:
#'   \code{ko} and a numeric column with the same name as the source. Missing
#'   source effects should be encoded as \code{NA}.
#' @param pname Name of the gene program node, such as \code{"PC1"}.
#' @param pgenes Character vector of gene program member genes.
#' @param ncores The number of cores for parallel computation.
#' @param alpha Significance level threshold (F-statistic p-value) used to
#'   determine whether a portability model is valid.
#' @return A list containing `B_matrix` and `portability_models`.
#'   `portability_models` contains a model for each combination of sources,
#'   adjusted R-squared values for model selection, and fallback values.
#' @export
fit_expression_model_with_gene_program <- function(
  Y, group, graph, source_effects, pname, pgenes, ncores, alpha = 0.1
) {
  .check_perturbation_effect_inputs(Y = Y, group = group)
  .check_ncores(ncores)
  .check_source_effects(source_effects)
  .check_alpha(alpha)

  nodes <- colnames(Y)
  .check_igraph(graph, nodes = nodes)
  .check_gene_program_nodes(pname = pname, pgenes = pgenes, nodes = nodes)

  # Stage 1: fit the target-dataset expression model.
  target_B_matrix <- fit_expression_model(
    Y = Y,
    group = group,
    graph = graph,
    method = 'lm',
    ncores = ncores
  )

  # Stage 2: train portability regressions for the gene program node.
  portability_models <- .fit_gene_program_portability_models(
    Y = Y,
    group = group,
    source_effects = source_effects,
    pname = pname,
    pgenes = pgenes,
    alpha = alpha
  )

  return(list(
    B_matrix = target_B_matrix,
    portability_models = portability_models
  ))
}


# Fit portability regressions for predicting target gene program node effects.
.fit_gene_program_portability_models <- function(
  Y, group, source_effects, pname, pgenes, alpha
) {
  nodes <- colnames(Y)
  training_kos <- setdiff(unique(group), 'WT')
  wt_means <- colMeans(Y[group == 'WT', , drop = FALSE])

  # Build the regression response from target data.
  # Portability response: observed target-dataset effect on the gene program
  # node only.  Source effects are the model inputs for this single response.
  target_program_effects <- dplyr::tibble(
    ko = training_kos,
    target_delta = vapply(
      training_kos,
      function(ko) {
        mean(Y[group == ko, pname]) - wt_means[pname]
      },
      FUN.VALUE = numeric(1)
    )
  )

  # Store fallback deltas for prediction when no valid source model is available.
  # Fallbacks are empirical target effects from the training KOs.
  training_ko_delta_matrix <- vapply(
    training_kos,
    function(ko) {
      colMeans(Y[group == ko, , drop = FALSE]) - wt_means
    },
    FUN.VALUE = setNames(numeric(length(nodes)), nodes)
  )
  mean_training_deltas <- rowMeans(training_ko_delta_matrix, na.rm = TRUE)

  source_names <- names(source_effects)
  portability_models <- list(
    fallback_pname_delta = mean(target_program_effects$target_delta, na.rm = TRUE),
    fallback_pgenes_deltas = mean_training_deltas[pgenes]
  )

  # Join the response and source effects into one training table keyed by KO.
  portability_training_df <- Reduce(
    function(x, y) dplyr::left_join(x, y, by = 'ko'),
    c(list(target_program_effects), source_effects)
  )

  source_combos <- unlist(
    lapply(seq_along(source_names), function(k) {
      utils::combn(source_names, k, simplify = FALSE)
    }),
    recursive = FALSE
  )

  # Train one model for each non-empty source combination. Prediction later uses
  # the model matching the sources available for a test KO.
  for (source_combo in source_combos) {
    model_name <- paste(sort(source_combo), collapse = "_and_")

    # Use only KOs with complete source-effect values for this source combination.
    complete_training_df <- portability_training_df[
      stats::complete.cases(portability_training_df[, source_combo, drop = FALSE]),
      ,
      drop = FALSE
    ]

    if (nrow(complete_training_df) > length(source_combo)) {
      portability_formula <- stats::as.formula(
        paste('target_delta ~', paste(source_combo, collapse = ' + '))
      )
      portability_fit <- stats::lm(portability_formula, data = complete_training_df)

      fit_summary <- summary(portability_fit)
      fstat <- fit_summary$fstatistic
      p_value <- stats::pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE)
      is_significant <- !is.na(p_value) && p_value < alpha
      adj_r_squared <- fit_summary$adj.r.squared

      # Linear portability models can get the direction right but shrink or
      # inflate the magnitude.  Store a scale correction so prediction can match
      # the observed target-program-node variance in the training KOs.
      observed_scale <- stats::sd(complete_training_df$target_delta)
      predicted_scale <- stats::sd(stats::predict(
        portability_fit,
        newdata = complete_training_df
      ))
      scale_factor <- ifelse(
        predicted_scale > 1e-6,
        observed_scale / predicted_scale,
        1.0
      )

      portability_models[[model_name]] <- list(
        model = portability_fit,
        scale_factor = scale_factor,
        is_significant = is_significant,
        adj_r_squared = adj_r_squared
      )
    } else {
      portability_models[[model_name]] <- list(
        model = NULL,
        scale_factor = 1.0,
        is_significant = FALSE,
        adj_r_squared = NA_real_
      )
    }
  }

  return(portability_models)
}


#' Predict Perturbation Effect With Gene Program
#'
#' Predicts the perturbation effect by matching test KOs to the best-available
#' pre-trained portability model and propagating the effect.
#'
#' @param B A numeric matrix of regulatory coefficients.
#' @param portability_models A named list of trained models and fallback values.
#' @param source_effects A named list of source-effect data frames in the same
#'   format used by \code{\link{fit_expression_model_with_gene_program}}.
#' @param ko_expressions A named vector of expression levels for the test KOs.
#' @param wt_expressions A named vector of WT expressions for the target dataset.
#' @param pname Name of the gene program node.
#' @param pgenes Character vector of gene program member genes.
#' @param scale_pname Logical. If TRUE, apply the stored variance calibration to
#'   the predicted gene program node effect.
#' @param model_selection Portability model selection strategy. \code{"exact"}
#'   preserves the original behavior and uses only the model matching all
#'   available sources for a KO. \code{"best_adj_r_squared"} chooses the
#'   significant available-source subset model with highest adjusted R-squared.
#' @return A numeric matrix of predicted delta values for the test KOs.
#' @export
predict_perturbation_effect_with_gene_program <- function(
  B, portability_models, source_effects, ko_expressions, wt_expressions,
  pname, pgenes, scale_pname = TRUE,
  model_selection = c('exact', 'best_adj_r_squared')
) {
  # Validate the target expression vectors, fitted GRN model, source effects, and
  # gene program definition before choosing any portability model.
  .check_expression_vector(wt_expressions, arg = 'wt_expressions')
  .check_expression_vector(ko_expressions, arg = 'ko_expressions')

  genes <- names(wt_expressions)
  ko_genes <- names(ko_expressions)

  .check_expression_model_matrix(B, genes = genes)
  .check_gene_program_nodes(pname = pname, pgenes = pgenes, nodes = genes)
  if (!all(ko_genes %in% genes)) {
    stop(
      "'ko_expressions' names must be included in 'wt_expressions'.",
      call. = FALSE
    )
  }
  .check_source_effects(source_effects)
  if (!is.list(portability_models)) {
    stop("'portability_models' must be a list.", call. = FALSE)
  }
  if (
    !is.logical(scale_pname) || length(scale_pname) != 1L || is.na(scale_pname)
  ) {
    stop("'scale_pname' must be TRUE or FALSE.", call. = FALSE)
  }
  model_selection <- match.arg(model_selection)

  B_propagator <- t(B[genes, , drop = FALSE])
  adj_matrix_functional <- t(B_propagator != 0)
  graph_for_distances <- igraph::graph_from_adjacency_matrix(
    adj_matrix_functional,
    mode = 'directed'
  )
  all_source_names <- names(source_effects)

  # Validate fallback values learned during fitting.  These are required when a
  # KO lacks source effects or its matching portability model is not significant.
  fallback_pname_delta <- portability_models[['fallback_pname_delta']]
  if (
    !is.numeric(fallback_pname_delta) || length(fallback_pname_delta) != 1L ||
      !is.finite(fallback_pname_delta)
  ) {
    stop(
      "Corrupted 'portability_models': 'fallback_pname_delta' is missing.",
      call. = FALSE
    )
  }

  fallback_pgenes_deltas <- portability_models[['fallback_pgenes_deltas']]
  if (
    !is.numeric(fallback_pgenes_deltas) ||
      is.null(names(fallback_pgenes_deltas))
  ) {
    stop(
      "Corrupted 'portability_models': 'fallback_pgenes_deltas' is missing.",
      call. = FALSE
    )
  }
  if (!setequal(pgenes, names(fallback_pgenes_deltas))) {
    stop("'pgenes' do not match 'fallback_pgenes_deltas'.", call. = FALSE)
  }

  pred_delta_list <- lapply(ko_genes, function(ko_gene) {
    use_fallback <- FALSE
    predicted_pname_delta <- NA_real_

    # Find which source datasets contain this KO effect.  The model name is the
    # sorted source combination, matching names created during fitting.
    available_sources <- all_source_names[vapply(
      all_source_names,
      function(source_name) {
        source_effect <- source_effects[[source_name]]
        source_value <- source_effect[[source_name]][source_effect$ko == ko_gene]
        length(source_value) == 1L && !is.na(source_value)
      },
      FUN.VALUE = logical(1)
    )]

    selected_sources <- character(0)
    model_entry <- NULL

    if (length(available_sources) == 0) {
      use_fallback <- TRUE
    } else if (model_selection == 'exact') {
      model_name <- paste(sort(available_sources), collapse = "_and_")
      model_entry <- portability_models[[model_name]]
      selected_sources <- available_sources

      # Use a portability model only if it exists and passed the fit-time
      # F-test threshold.  Otherwise fall back to target training averages.
      if (is.null(model_entry) || !isTRUE(model_entry$is_significant)) {
        use_fallback <- TRUE
      }
    } else {
      # Consider every trained model whose source set is available for this KO.
      source_subsets <- unlist(
        lapply(seq_along(available_sources), function(k) {
          utils::combn(available_sources, k, simplify = FALSE)
        }),
        recursive = FALSE
      )
      candidate_model_names <- vapply(
        source_subsets,
        function(source_subset) {
          paste(sort(source_subset), collapse = "_and_")
        },
        FUN.VALUE = character(1)
      )
      candidate_models <- portability_models[candidate_model_names]

      valid_candidates <- vapply(
        candidate_models,
        function(candidate_model) {
          !is.null(candidate_model) && !is.null(candidate_model$model) &&
            isTRUE(candidate_model$is_significant)
        },
        FUN.VALUE = logical(1)
      )

      if (!any(valid_candidates)) {
        use_fallback <- TRUE
      } else {
        valid_models <- candidate_models[valid_candidates]
        valid_source_subsets <- source_subsets[valid_candidates]
        adj_r_squared <- vapply(
          valid_models,
          function(valid_model) {
            if (
              !is.numeric(valid_model$adj_r_squared) ||
                length(valid_model$adj_r_squared) != 1L ||
                is.na(valid_model$adj_r_squared)
            ) {
              stop(
                "Portability models are missing 'adj_r_squared'; refit with ",
                "fit_expression_model_with_gene_program().",
                call. = FALSE
              )
            }
            valid_model$adj_r_squared
          },
          FUN.VALUE = numeric(1)
        )
        selected_idx <- which.max(adj_r_squared)
        model_entry <- valid_models[[selected_idx]]
        selected_sources <- valid_source_subsets[[selected_idx]]
      }
    }

    if (!use_fallback) {
      # Build the single-row source-effect table for this KO, then predict the
      # target gene program node delta.
      ko_data_list <- lapply(selected_sources, function(source_name) {
        dplyr::filter(source_effects[[source_name]], ko == ko_gene)
      })
      ko_data <- Reduce(
        function(x, y) dplyr::full_join(x, y, by = 'ko'),
        ko_data_list
      )

      raw_pred <- stats::predict(model_entry$model, newdata = ko_data)

      predicted_pname_delta <- if (scale_pname) {
        raw_pred * model_entry$scale_factor
      } else {
        raw_pred
      }
    }

    if (use_fallback) {
      predicted_pname_delta <- fallback_pname_delta
    }

    # Known effects before propagation: the measured KO-gene delta in the target
    # dataset and the predicted or fallback effect on the gene program node.
    known_deltas <- numeric(0)
    known_deltas[ko_gene] <- ko_expressions[ko_gene] - wt_expressions[ko_gene]
    known_deltas[pname] <- predicted_pname_delta

    # When the program-node effect is a fallback, also anchor member genes to
    # their target training averages rather than imputing them only from pname.
    if (use_fallback) {
      pgenes_to_set <- setdiff(pgenes, ko_gene)
      if (length(pgenes_to_set) > 0) {
        known_deltas[pgenes_to_set] <- fallback_pgenes_deltas[pgenes_to_set]
      }
    }

    distances <- igraph::distances(
      graph_for_distances,
      v = names(known_deltas),
      mode = 'out'
    )
    unreachable_genes <- names(which(colSums(is.finite(distances)) == 0L))
    known_deltas[unreachable_genes] <- 0

    # Propagate the known deltas through the fitted target GRN by solving for the
    # steady-state deltas of all remaining genes.
    imputed_deltas <- .impute_deltas(
      B_propagator = B_propagator,
      known_deltas = known_deltas
    )
    if (length(imputed_deltas) > 0) {
      known_deltas[names(imputed_deltas)] <- imputed_deltas
    }

    if (!setequal(names(known_deltas), genes)) {
      stop("Prediction did not return a value for all genes.", call. = FALSE)
    }
    if (anyNA(known_deltas)) {
      stop("Prediction resulted in missing values.", call. = FALSE)
    }

    return(known_deltas[genes])
  })

  # Assemble one predicted delta vector per KO.
  pred_matrix <- do.call(rbind, pred_delta_list)
  rownames(pred_matrix) <- ko_genes

  return(pred_matrix)
}
