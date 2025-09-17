#' Fit Expression Model with a Gene Program
#'
#' This function prepares the two key components needed for the sophisticated
#' hybrid prediction model. It leverages the standard \code{\link{fit_expression_model}}
#' to learn the main linear propagator (B matrix) after applying specific graph
#' modifications, and separately learns the robust, univariate "switch"
#' coefficients for the hub's triggers.
#'
#' @param Y A numeric matrix of expression data (cells x genes), including the hub node.
#' @param group A character or factor vector indicating the group for each cell in Y.
#' @param graph An igraph object representing the full regulatory structure.
#' @param pname A character string for the name of the hub node (e.g., "PC1").
#' @param pgenes A character vector of the signature genes regulated by the hub.
#' @param ncores The number of cores for parallel computation.
#' @param method The regression method for the main B matrix. One of 'lm', 'lasso', or 'ridge'.
#' @param pgene_parent_set A character string specifying the allowed set of
#'   regulators for signature genes. One of 'all', 'non_signature', or 'hub_only'.
#' @return A list with two components:
#'   \item{B_matrix}{The main regulatory matrix (sources x targets).}
#'   \item{switch_coefficients}{A list where each element is the named coefficient
#'     vector from the univariate regression of the hub against one of its parents.}
#' @export
fit_expression_model_with_gp <- function(Y, group, graph, pname, pgenes, ncores,
                                         method = c('ridge', 'lasso', 'lm'),
                                         pgene_parent_set = c('all', 'non_signature', 'hub_only')) {

  # --- 1. Parameter Validation ---
  all_nodes <- colnames(Y)
  method <- match.arg(method)
  pgene_parent_set <- match.arg(pgene_parent_set)

  stopifnot(inherits(graph, "igraph"))
  stopifnot(setequal(all_nodes, igraph::V(graph)$name))
  stopifnot(!igraph::any_loop(graph))
  stopifnot(all(c(pname, pgenes) %in% all_nodes))

  message("--- Preparing components for hybrid model ---")

  # --- 2. Learn the Univariate "Switch" Coefficients ---
  message("Step 1/2: Learning univariate 'switch' coefficients for the hub...")
  hub_parents <- igraph::neighbors(graph, pname, mode = 'in')$name

  if (length(hub_parents) > 0) {
    switch_coefficients <- pbmcapply::pbmclapply(hub_parents, function(parent) {
      samples <- which(group != parent & group != pname) # A conservative sample set
      fit <- stats::lm(Y[samples, pname] ~ Y[samples, parent])
      return(stats::coef(fit))
    }, mc.cores = ncores)
    names(switch_coefficients) <- hub_parents
  } else {
    # If the hub has no parents, the model is misspecified.
    stop(paste("Hub node '", pname, "' has no parents in the graph. Consider removing it before analysis."))
  }

  # --- 3. Learn the Main B Matrix (The "Propagator") ---
  message(paste0("Step 2/2: Learning main B matrix with pgene_parent_set = '", pgene_parent_set, "'..."))

  # Create a modified graph based on the pgene_parent_set rule
  graph_for_B <- graph
  if (pgene_parent_set != "all") {

    # Identify the set of predictors to disallow
    if (pgene_parent_set == "non_signature") {
      # Disallow other signature genes from being parents
      disallowed_parents <- pgenes
    } else { # hub_only
      # Disallow everything except the hub from being a parent
      disallowed_parents <- setdiff(all_nodes, pname)
    }

    adj_matrix <- igraph::as_adjacency_matrix(graph, sparse = FALSE)
    # Set the disallowed links to 0
    adj_matrix[disallowed_parents, pgenes] <- 0
    # Convert back to an igraph object
    graph_for_B <- igraph::graph_from_adjacency_matrix(adj_matrix, mode = 'directed')
  }

  # Call the existing, robust function to fit the B matrix
  B_matrix <- fit_expression_model(
    Y = Y,
    group = group,
    graph = graph_for_B,
    ncores = ncores,
    method = method
  )

  message("--- Component training complete. ---")

  # --- 4. Return Both Components ---
  return(list(
    B_matrix = B_matrix,
    switch_coefficients = switch_coefficients
  ))
}


#' Predict Perturbation Effect with a Gene Program Switch
#'
#' This function implements a two-stage prediction method in presence of a gene program switch.
#' It partitions the network and solves for components in a logical order to
#' correctly handle non-linear "switch" effects at a central hub (pname).
#'
#' @param B A numeric matrix of regulatory coefficients (sources x targets),
#'   including the hub node. This is the output of \code{\link{fit_expression_model_with_gp}}.
#' @param switch_coefs A named list of univariate regression coefficients for the
#'   hub's parents. The function assumes that all `pgenes` are parents of the
#'   hub and thus present as names in this list.
#' @param ko_expressions A named vector of expression levels for perturbed genes.
#' @param wt_expressions A named vector of WT expressions (including the hub node).
#' @param pname The name of the hub node (e.g., "PC1").
#' @param pgenes A character vector of the signature genes.
#' @param pdirection An integer (-1 or 1) specifying the consistent direction
#'   of the switch effect on the hub. Set as 1 if pname increases under perturbation,
#'   -1 if pname decreases under perturbation.
#' @param max_dist Maximum distance for effect propagation. Default is `Inf`.
#'   Distances are calculated on the functional graph derived from B.
#' @param space The modeling space. delta: propagate in expression change space. absolute: propagate in expression space.
#' @param ncores The number of cores for parallel computation.
#' @return A numeric matrix of predicted delta values (KOs x genes, including pname).
#' @export
predict_standard_effect_with_gp <- function(
    B, switch_coefs, ko_expressions, wt_expressions, pname, pgenes, pdirection,
    max_dist = Inf, space = c("delta", "absolute"), ncores = 1
) {
  # --- 1. Input Validation and Setup ---
  space <- match.arg(space)
  genes <- names(wt_expressions)
  ko_genes <- names(ko_expressions)
  stopifnot(
    is.matrix(B), is.numeric(wt_expressions), is.numeric(ko_expressions),
    identical(colnames(B), genes),
    identical(rownames(B), c('Intercept', genes)),
    all(c(pname, pgenes, ko_genes) %in% genes),
    pdirection %in% c(-1, 1),
    max_dist >= 2,
    all(pgenes %in% names(switch_coefs))
  )

  non_signature_genes <- setdiff(genes, pgenes)

  # --- 2. Prepare B Matrix and Functional Graph ---
  B_propagator <- t(B[genes, , drop = FALSE])
  adj_matrix_functional <- t(B_propagator != 0)
  graph_for_distances <- igraph::graph_from_adjacency_matrix(adj_matrix_functional, mode = 'directed')
  all_pname_parents <- igraph::neighbors(graph_for_distances, pname, mode = 'in')$name
  B_propagator_nonsig <- B_propagator[non_signature_genes, non_signature_genes, drop = FALSE]

  all_distances <- igraph::distances(graph_for_distances, v = ko_genes, mode = 'out')

  # --- 3. Main Prediction Loop ---
  pred_delta_list <- pbmclapply::pbmclapply(ko_genes, function(ko_gene) {

    distances <- all_distances[ko_gene, ]
    unchanged_genes <- names(which(is.infinite(distances) | distances > max_dist))

    if (ko_gene %in% pgenes) {
      # --- CASE A: KO is a SIGNATURE gene ---
      if (space == 'delta') {
        known_deltas <- c()
        known_deltas[ko_gene] <- ko_expressions[ko_gene] - wt_expressions[ko_gene]
        known_deltas[unchanged_genes] <- 0

        stopifnot(!(pname %in% unchanged_genes))
        known_deltas[pname] <- -pdirection * abs(switch_coefs[[ko_gene]][2]) * known_deltas[ko_gene]

        imputed_deltas <- .impute_deltas(B_propagator = B_propagator, known_deltas = known_deltas)
        if (length(imputed_deltas) > 0) known_deltas[names(imputed_deltas)] <- imputed_deltas

        stopifnot(length(known_deltas) == length(genes), setequal(names(known_deltas), genes))
        return(known_deltas[genes])
      } else { # space == "absolute"
        pred_ko_abs <- c()
        pred_ko_abs[ko_gene] <- ko_expressions[ko_gene]
        pred_ko_abs[unchanged_genes] <- wt_expressions[unchanged_genes]

        stopifnot(!(pname %in% unchanged_genes))
        delta_pname <- -pdirection * abs(switch_coefs[[ko_gene]][2]) * (ko_expressions[ko_gene] - wt_expressions[ko_gene])
        pred_ko_abs[pname] <- wt_expressions[pname] + delta_pname

        imputed_abs <- .impute_absolute(B_full = t(B), known_expressions = pred_ko_abs)
        if (length(imputed_abs) > 0) pred_ko_abs[names(imputed_abs)] <- imputed_abs

        stopifnot(length(pred_ko_abs) == length(genes), setequal(names(pred_ko_abs), genes))
        return(pred_ko_abs[genes] - wt_expressions)
      }
    } else {
      # --- CASE B: KO is a NON-signature gene ---
      if (space == 'delta') {
        known_deltas <- c()
        known_deltas[ko_gene] <- ko_expressions[ko_gene] - wt_expressions[ko_gene]
        known_deltas[unchanged_genes] <- 0

        sub_known_deltas <- known_deltas[intersect(names(known_deltas), non_signature_genes)]
        imputed_deltas_nonsig <- .impute_deltas(B_propagator = B_propagator_nonsig, known_deltas = sub_known_deltas)
        if (length(imputed_deltas_nonsig) > 0) known_deltas[names(imputed_deltas_nonsig)] <- imputed_deltas_nonsig

        pname_parents_nonsig <- intersect(all_pname_parents, non_signature_genes)
        parent_deltas <- known_deltas[pname_parents_nonsig]
        if (length(parent_deltas) > 0 && any(parent_deltas != 0, na.rm = TRUE)) {
          deltas_pname <- sapply(pname_parents_nonsig, function(p) abs(switch_coefs[[p]][2]) * abs(parent_deltas[p]))
          strongest_source <- names(which.max(deltas_pname))
          delta_pname <- -pdirection * abs(switch_coefs[[strongest_source]][2]) * known_deltas[strongest_source]
        } else {
          delta_pname <- 0
        }
        known_deltas[pname] <- delta_pname

        imputed_deltas_sig <- .impute_deltas(B_propagator = B_propagator, known_deltas = known_deltas)
        if (length(imputed_deltas_sig) > 0) known_deltas[names(imputed_deltas_sig)] <- imputed_deltas_sig

        stopifnot(length(known_deltas) == length(genes), setequal(names(known_deltas), genes))
        return(known_deltas[genes])
      } else { # space == "absolute"
        pred_ko_abs <- c()
        pred_ko_abs[ko_gene] <- ko_expressions[ko_gene]
        pred_ko_abs[unchanged_genes] <- wt_expressions[unchanged_genes]

        sub_pred_ko_abs <- pred_ko_abs[intersect(names(pred_ko_abs), non_signature_genes)]
        imputed_abs_nonsig <- .impute_absolute(
          B_full = t(B)[non_signature_genes, c('Intercept', non_signature_genes)],
          known_expressions = sub_pred_ko_abs
        )
        if (length(imputed_abs_nonsig) > 0) pred_ko_abs[names(imputed_abs_nonsig)] <- imputed_abs_nonsig

        pname_parents_nonsig <- intersect(all_pname_parents, non_signature_genes)
        parent_deltas <- pred_ko_abs[pname_parents_nonsig] - wt_expressions[pname_parents_nonsig]
        if (length(parent_deltas) > 0 && any(parent_deltas != 0, na.rm = TRUE)) {
          deltas_pname <- sapply(pname_parents_nonsig, function(p) abs(switch_coefs[[p]][2]) * abs(parent_deltas[p]))
          strongest_source <- names(which.max(deltas_pname))
          delta_strongest_source <- pred_ko_abs[strongest_source] - wt_expressions[strongest_source]
          delta_pname <- -pdirection * abs(switch_coefs[[strongest_source]][2]) * delta_strongest_source
        } else {
          delta_pname <- 0
        }
        pred_ko_abs[pname] <- wt_expressions[pname] + delta_pname

        imputed_abs <- .impute_absolute(B_full = t(B), known_expressions = pred_ko_abs)
        if (length(imputed_abs) > 0) pred_ko_abs[names(imputed_abs)] <- imputed_abs

        stopifnot(length(pred_ko_abs) == length(genes), setequal(names(pred_ko_abs), genes))
        return(pred_ko_abs[genes] - wt_expressions)
      }
    }
  }, mc.cores = ncores)

  # --- 4. Assemble and Return Final Matrix ---
  pred_matrix <- do.call(rbind, pred_delta_list)
  rownames(pred_matrix) <- ko_genes

  return(pred_matrix)
}
