#' Fit Expression Model and Learn Average SCOPE Effect
#'
#' Learns the linear GRN model (B matrix) and calculates the average observed
#' effect on the SCOPE signature from the provided training data.
#'
#' @param Y A numeric matrix of expression data (cells x genes), including the hub node.
#' @param group A character or factor vector indicating the group for each cell in Y.
#' @param graph An igraph object representing the full regulatory structure.
#' @param pname A character string for the name of the hub node (e.g., "PC1").
#' @param pgenes A character vector of the signature genes regulated by the hub.
#' @param de_true A ground truth matrix of DE values (KOs x genes) for the
#'   training set, used to calculate the average SCOPE effect.
#' @param stat A long-format data frame of statistics for the training set,
#'   used to identify KOs that trigger the hub.
#' @param ncores The number of cores for parallel computation.
#' @param method The regression method for the main B matrix. One of 'lm', 'lasso', or 'ridge'.
#' @param pname_to_pgenes_only A logical to constrain hub effects to signature genes.
#' @return A list with two components:
#'   \item{B_matrix}{The main regulatory matrix (sources x targets).}
#'   \item{scope_mean_effect}{A named vector of the average DE for the hub and
#'     signature genes, calculated from training KOs that trigger the hub.}
#' @export
fit_expression_model_with_gp <- function(
    Y, group, graph, pname, pgenes, de_true, stat, ncores,
    method = c('lm', 'ridge', 'lasso'),
    pname_to_pgenes_only = TRUE
) {

  # --- 1. Parameter Validation ---
  all_nodes <- colnames(Y)
  method <- match.arg(method)

  stopifnot(
    all(c(pname, pgenes) %in% colnames(de_true)),
    inherits(graph, "igraph"),
    setequal(all_nodes, igraph::V(graph)$name),
    !igraph::any_loop(graph),
    all(c(pname, pgenes) %in% all_nodes)
  )

  message("--- Preparing components for hybrid model ---")

  # --- 2. Learn the Average SCOPE Effect from Training Data ---
  message("Step 1/2: Learning average SCOPE effect from training data...")
  kos_triggering_pname <- stat %>%
    dplyr::filter(gene == pname, adj_pv < 0.05) %>%
    dplyr::pull(ko) %>%
    unique()

  scope_signature <- c(pname, pgenes)

  if (length(kos_triggering_pname) > 0) {
    scope_mean_effect <- colMeans(
      de_true[intersect(rownames(de_true), kos_triggering_pname), scope_signature, drop = FALSE],
      na.rm = TRUE
    )
  } else {
    warning("No training KOs found that significantly trigger the hub. Average SCOPE effect will be zero.")
    scope_mean_effect <- setNames(rep(0, length(scope_signature)), scope_signature)
  }

  # --- 3. Learn the Main B Matrix (The "Propagator") ---
  message("Step 2/2: Learning main B matrix...")
  graph_for_B <- graph
  if (pname_to_pgenes_only) {
    adj_matrix_for_B <- igraph::as_adjacency_matrix(graph_for_B, sparse = FALSE)
    non_signature_genes <- setdiff(all_nodes, c(pname, pgenes))
    adj_matrix_for_B[pname, non_signature_genes] <- 0
    graph_for_B <- igraph::graph_from_adjacency_matrix(adj_matrix_for_B, mode = 'directed')
  }

  B_matrix <- fit_expression_model(
    Y = Y, group = group, graph = graph_for_B,
    ncores = ncores, method = method
  )

  message("--- Component training complete. ---")

  # --- 4. Return Both Components ---
  return(list(
    B_matrix = B_matrix,
    scope_mean_effect = scope_mean_effect
  ))
}


#' Predict Perturbation Effect Using a Data-Driven SCOPE Imputation
#'
#' Predicts KO effects using a hybrid model. If a KO is predicted to affect the
#' SCOPE signature, the effect on the hub and signature genes is imputed using a
#' pre-calculated average effect. Effects on non-signature genes are then
#' calculated using a linear solver.
#'
#' @param B A numeric matrix of regulatory coefficients (sources x targets).
#' @param scope_mean_effect A named vector of the average training DE for the hub
#'   and signature genes.
#' @param ko_expressions A named vector of expression levels for perturbed genes.
#' @param wt_expressions A named vector of WT expressions.
#' @param pname The name of the hub node (e.g., "PC1").
#' @param pgenes A character vector of the signature genes.
#' @param max_dist Maximum distance for effect propagation.
#' @param ncores Number of cores for parallel computation.
#' @return A numeric matrix of predicted delta values (KOs x genes).
#' @export
predict_standard_effect_with_gp <- function(
    B, scope_mean_effect, ko_expressions, wt_expressions, pname, pgenes,
    max_dist = Inf, ncores = 1
) {
  # --- 1. Input Validation and Setup ---
  genes <- names(wt_expressions)
  ko_genes <- names(ko_expressions)

  stopifnot(
    is.matrix(B), is.numeric(wt_expressions), is.numeric(ko_expressions),
    identical(colnames(B), genes),
    identical(rownames(B), c('Intercept', genes)),
    all(c(pname, pgenes, ko_genes) %in% genes),
    max_dist >= 1 # A distance of 1 is the minimum meaningful value
  )

  # --- 2. Prepare B Matrix and Functional Graph ---
  B_propagator <- t(B[genes, , drop = FALSE])
  adj_matrix_functional <- t(B_propagator != 0)
  graph_for_distances <- igraph::graph_from_adjacency_matrix(adj_matrix_functional, mode = 'directed')
  all_distances <- igraph::distances(graph_for_distances, v = ko_genes, mode = 'out')

  # --- 3. Main Prediction Loop ---
  pred_delta_list <- pbmcapply::pbmclapply(ko_genes, function(ko_gene) {

    distances <- all_distances[ko_gene, ]
    unchanged_genes <- names(which(is.infinite(distances) | distances > max_dist))

    known_deltas <- c()
    known_deltas[ko_gene] <- ko_expressions[ko_gene] - wt_expressions[ko_gene]
    known_deltas[unchanged_genes] <- 0

    is_pname_triggered <- !(pname %in% unchanged_genes)

    if (is_pname_triggered) {
      scope_genes_to_impute <- setdiff(names(scope_mean_effect), ko_gene)
      known_deltas[scope_genes_to_impute] <- scope_mean_effect[scope_genes_to_impute]
    }

    imputed_deltas <- .impute_deltas(B_propagator = B_propagator, known_deltas = known_deltas)
    if (length(imputed_deltas) > 0) known_deltas[names(imputed_deltas)] <- imputed_deltas

    stopifnot(
      length(known_deltas) == length(genes),
      setequal(names(known_deltas), genes),
      !any(is.na(known_deltas))
    )

    return(known_deltas[genes])

  }, mc.cores = ncores)

  # --- 4. Assemble and Return Final Matrix ---
  pred_matrix <- do.call(rbind, pred_delta_list)
  rownames(pred_matrix) <- ko_genes

  return(pred_matrix)
}
