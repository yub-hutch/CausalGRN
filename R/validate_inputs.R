# Input validation helpers
#
# Shared checks:
# - .check_ncores(): check positive core count.
# - .check_group(): check named character group labels and optional cell counts.
# - .check_count_matrix(): check dense, finite, non-negative count matrix.
# - .check_expression_matrix(): check finite numeric expression matrix.
# - .check_adjacency_matrix(): check optional adjacency matrix against nodes.
#
# Function-specific checks:
# - .check_expression_simulator_params(): check simulator scalar parameters.
# - .check_perturbation_effect_inputs(): check Y and group consistency.
# - .check_perturbation_effect_params(): check perturbation-effect parameters.
# - .check_skeleton_data(): check count, Y, and G consistency.
# - .check_skeleton_params(): check skeleton scalar tuning parameters.
# - .check_causalgrn_data(): check graph and perturbation statistics.
# - .check_causalgrn_params(): check causal orientation scalar parameters.


# Check core count input.
.check_ncores <- function(ncores) {
  if (
    length(ncores) != 1L || !is.numeric(ncores) ||
      is.na(ncores) || !is.finite(ncores) ||
      ncores < 1 || ncores != round(ncores)
  ) {
    stop("'ncores' must be a single positive integer.", call. = FALSE)
  }

  invisible(TRUE)
}


# Check group labels.
.check_group <- function(group, min_cells = NULL) {
  if (!is.character(group) || !is.null(dim(group))) {
    stop("'group' must be a character vector.", call. = FALSE)
  }
  if (is.null(names(group))) {
    stop("'group' must be named.", call. = FALSE)
  }
  if (anyNA(names(group)) || any(names(group) == "")) {
    stop("'group' names must be non-missing.", call. = FALSE)
  }
  if (anyDuplicated(names(group))) {
    stop("'group' names must be unique.", call. = FALSE)
  }

  if (any(is.na(group) | group == "")) {
    stop("'group' must not contain missing or empty labels.", call. = FALSE)
  }
  if (!"WT" %in% group) {
    stop("'group' must contain wild-type cells labeled as 'WT'.", call. = FALSE)
  }
  if (!is.null(min_cells)) {
    if (
      length(min_cells) != 1L || !is.numeric(min_cells) ||
        !is.finite(min_cells) || min_cells < 10 || min_cells != round(min_cells)
    ) {
      stop("'min_cells' must be a single integer >= 10.", call. = FALSE)
    }
    if (any(table(group) < min_cells)) {
      stop("Each group must contain at least ", min_cells, " cells.", call. = FALSE)
    }
  }

  invisible(TRUE)
}


# Check that count data are a dense, finite, non-negative, integer-like matrix.
.check_count_matrix <- function(count) {
  if (!is.matrix(count) || !is.numeric(count)) {
    stop("'count' must be a numeric matrix.", call. = FALSE)
  }
  if (is.null(rownames(count)) || is.null(colnames(count))) {
    stop("'count' must have row and column names.", call. = FALSE)
  }
  if (
    anyNA(rownames(count)) || anyNA(colnames(count)) ||
      any(rownames(count) == "") || any(colnames(count) == "")
  ) {
    stop("'count' row and column names must be non-missing.", call. = FALSE)
  }
  if (anyDuplicated(rownames(count)) || anyDuplicated(colnames(count))) {
    stop("'count' row and column names must be unique.", call. = FALSE)
  }

  tol <- sqrt(.Machine$double.eps)
  for (j in seq_len(ncol(count))) {
    x <- count[, j]
    if (any(!is.finite(x))) {
      stop("'count' must contain only finite values.", call. = FALSE)
    }
    if (any(x < 0)) {
      stop("'count' must contain only non-negative values.", call. = FALSE)
    }
    if (any(abs(x - round(x)) > tol)) {
      stop("'count' must contain integer-like values.", call. = FALSE)
    }
  }

  invisible(TRUE)
}


# Check that normalized expression data are a finite numeric matrix.
.check_expression_matrix <- function(Y, arg = "Y") {
  arg <- paste0("'", arg, "'")
  if (!is.matrix(Y) || !is.numeric(Y)) {
    stop(arg, " must be a numeric matrix.", call. = FALSE)
  }
  if (is.null(rownames(Y)) || is.null(colnames(Y))) {
    stop(arg, " must have row and column names.", call. = FALSE)
  }
  if (
    anyNA(rownames(Y)) || anyNA(colnames(Y)) ||
      any(rownames(Y) == "") || any(colnames(Y) == "")
  ) {
    stop(arg, " row and column names must be non-missing.", call. = FALSE)
  }
  if (anyDuplicated(rownames(Y)) || anyDuplicated(colnames(Y))) {
    stop(arg, " row and column names must be unique.", call. = FALSE)
  }

  for (j in seq_len(ncol(Y))) {
    if (any(!is.finite(Y[, j]))) {
      stop(arg, " must contain only finite values.", call. = FALSE)
    }
  }

  invisible(TRUE)
}


# Check an adjacency matrix against an exact node order.
.check_adjacency_matrix <- function(G, nodes) {
  if (is.null(G)) {
    return(invisible(TRUE))
  }

  if (!is.matrix(G) || nrow(G) != ncol(G)) {
    stop("'G' must be a square adjacency matrix.", call. = FALSE)
  }
  if (!(is.logical(G) || is.numeric(G))) {
    stop("'G' must be a logical or numeric adjacency matrix.", call. = FALSE)
  }
  if (is.null(rownames(G)) || is.null(colnames(G))) {
    stop("'G' must have row and column names.", call. = FALSE)
  }
  if (!identical(rownames(G), nodes) || !identical(colnames(G), nodes)) {
    stop("'G' row and column names must exactly match expected nodes.", call. = FALSE)
  }
  if (anyNA(G)) {
    stop("'G' must not contain missing values.", call. = FALSE)
  }

  G_logical <- G != 0
  if (any(diag(G_logical))) {
    stop("'G' must not contain self-loops.", call. = FALSE)
  }
  invisible(TRUE)
}


# Check scalar parameters for GRN-guided expression simulation.
.check_expression_simulator_params <- function(
    d, min_coef, max_coef, center_normal_sd, center_ko_eff, max_attempts
) {
  if (length(d) != 1L || !is.numeric(d) || !is.finite(d) || d <= 1) {
    stop("'d' must be a single finite number > 1.", call. = FALSE)
  }

  if (
    length(min_coef) != 1L || !is.numeric(min_coef) ||
      !is.finite(min_coef) || min_coef < 0
  ) {
    stop("'min_coef' must be a single finite non-negative number.", call. = FALSE)
  }

  if (
    length(max_coef) != 1L || !is.numeric(max_coef) ||
      !is.finite(max_coef) || max_coef < min_coef
  ) {
    stop("'max_coef' must be a single finite number >= 'min_coef'.", call. = FALSE)
  }

  if (
    length(center_normal_sd) != 1L ||
      !is.numeric(center_normal_sd) ||
      !is.finite(center_normal_sd) ||
      center_normal_sd <= 0
  ) {
    stop("'center_normal_sd' must be a single finite positive number.", call. = FALSE)
  }

  if (
    length(center_ko_eff) != 1L ||
      !is.numeric(center_ko_eff) ||
      !is.finite(center_ko_eff) ||
      center_ko_eff <= 0 ||
      center_ko_eff >= 1
  ) {
    stop("'center_ko_eff' must be a single finite number between 0 and 1.", call. = FALSE)
  }

  if (
    length(max_attempts) != 1L ||
      !is.numeric(max_attempts) ||
      !is.finite(max_attempts) ||
      max_attempts < 1 ||
      max_attempts != round(max_attempts)
  ) {
    stop("'max_attempts' must be a single positive integer.", call. = FALSE)
  }

  invisible(TRUE)
}


# Check inputs for perturbation-effect calculation.
.check_perturbation_effect_inputs <- function(Y, group) {
  if (length(group) != nrow(Y)) {
    stop("'group' must have one label per row of 'Y'.", call. = FALSE)
  }
  if (!identical(names(group), rownames(Y))) {
    stop("'group' names must be identical to 'Y' row names.", call. = FALSE)
  }

  kos <- setdiff(group, "WT")
  if (!length(kos)) {
    stop("'group' must contain at least one perturbation label.", call. = FALSE)
  }

  missing_kos <- setdiff(kos, colnames(Y))
  if (length(missing_kos)) {
    stop("Perturbed genes in 'group' must be columns of 'Y'.", call. = FALSE)
  }

  invisible(TRUE)
}


# Check scalar parameters for perturbation-effect calculation.
.check_perturbation_effect_params <- function(gene_block_size) {
  if (
    !is.null(gene_block_size) &&
      (
        length(gene_block_size) != 1L ||
          !is.numeric(gene_block_size) ||
          !is.finite(gene_block_size) ||
          gene_block_size < 1 ||
          gene_block_size != round(gene_block_size)
      )
  ) {
    stop("'gene_block_size' must be NULL or a single positive integer.", call. = FALSE)
  }

  invisible(TRUE)
}


# Check count, normalized expression, and optional graph inputs for skeleton inference.
.check_skeleton_data <- function(count, Y, G = NULL) {
  .check_count_matrix(count)
  .check_expression_matrix(Y)
  if (!identical(dim(count), dim(Y))) {
    stop("'count' and 'Y' must have identical dimensions.", call. = FALSE)
  }
  if (!identical(colnames(count), colnames(Y))) {
    stop("'count' and 'Y' must have identical column names.", call. = FALSE)
  }
  if (!identical(rownames(count), rownames(Y))) {
    stop("'count' and 'Y' must have identical row names.", call. = FALSE)
  }
  if (nrow(count) < 3L) {
    stop("'count' and 'Y' must contain at least 3 rows for skeleton inference.", call. = FALSE)
  }

  .check_adjacency_matrix(G, nodes = colnames(Y))

  invisible(TRUE)
}


# Check scalar tuning parameters for skeleton inference.
.check_skeleton_params <- function(
    alpha, min_abspcor, max_order, max_thr, min_n1, min_n2, sepset
) {
  if (
    length(alpha) != 1L || !is.numeric(alpha) ||
      !is.finite(alpha) || alpha <= 0 || alpha >= 1
  ) {
    stop("'alpha' must be a single finite number between 0 and 1.", call. = FALSE)
  }
  if (
    length(min_abspcor) != 1L || !is.numeric(min_abspcor) ||
      !is.finite(min_abspcor) || min_abspcor < 0 || min_abspcor > 1
  ) {
    stop("'min_abspcor' must be a single finite number between 0 and 1.", call. = FALSE)
  }
  if (
    length(max_order) != 1L || !is.numeric(max_order) ||
      !is.finite(max_order) || max_order != round(max_order) ||
      !(max_order %in% c(0, 1))
  ) {
    stop("'max_order' must be 0 or 1.", call. = FALSE)
  }
  if (
    length(max_thr) != 1L || !is.numeric(max_thr) ||
      !is.finite(max_thr) || max_thr < -1 || max_thr != round(max_thr)
  ) {
    stop("'max_thr' must be a single integer >= -1.", call. = FALSE)
  }
  if (
    length(min_n1) != 1L || !is.numeric(min_n1) ||
      !is.finite(min_n1) || min_n1 < 1 || min_n1 != round(min_n1)
  ) {
    stop("'min_n1' must be a single positive integer.", call. = FALSE)
  }
  if (
    length(min_n2) != 1L || !is.numeric(min_n2) ||
      !is.finite(min_n2) || min_n2 < 1 || min_n2 != round(min_n2)
  ) {
    stop("'min_n2' must be a single positive integer.", call. = FALSE)
  }
  if (!is.logical(sepset) || length(sepset) != 1L || is.na(sepset)) {
    stop("'sepset' must be TRUE or FALSE.", call. = FALSE)
  }

  invisible(TRUE)
}


# Check graph and perturbation statistics for causal GRN orientation.
.check_causalgrn_data <- function(graph, stat) {
  if (!inherits(graph, "igraph")) {
    stop("'graph' must be an igraph object.", call. = FALSE)
  }

  nodes <- igraph::V(graph)$name
  if (is.null(nodes) || anyNA(nodes) || any(nodes == "")) {
    stop("'graph' vertices must have non-missing names.", call. = FALSE)
  }
  if (anyDuplicated(nodes)) {
    stop("'graph' vertex names must be unique.", call. = FALSE)
  }

  required_cols <- c("ko", "gene", "adj_pv")
  if (!is.data.frame(stat) || !all(required_cols %in% colnames(stat))) {
    stop(
      "'stat' must be a data frame with columns 'ko', 'gene', and 'adj_pv'.",
      call. = FALSE
    )
  }
  if (
    anyNA(stat$ko) || anyNA(stat$gene) ||
      any(stat$ko == "") || any(stat$gene == "")
  ) {
    stop("'stat$ko' and 'stat$gene' must be non-missing names.", call. = FALSE)
  }
  if (!is.numeric(stat$adj_pv) || anyNA(stat$adj_pv)) {
    stop("'stat$adj_pv' must be numeric and non-missing.", call. = FALSE)
  }
  if (any(stat$adj_pv < 0 | stat$adj_pv > 1)) {
    stop("'stat$adj_pv' values must be between 0 and 1.", call. = FALSE)
  }
  if (anyDuplicated(paste(stat$ko, stat$gene, sep = "\r"))) {
    stop("'stat' must contain at most one row for each ko-gene pair.", call. = FALSE)
  }

  kos <- unique(stat$ko)
  genes <- unique(stat$gene)
  if (!all(kos %in% genes)) {
    stop(
      "All perturbed genes in 'stat$ko' must also appear in 'stat$gene'.",
      call. = FALSE
    )
  }
  if (nrow(stat) != length(kos) * length(genes)) {
    stop("'stat' must contain all ko-gene combinations exactly once.", call. = FALSE)
  }
  if (!setequal(genes, nodes)) {
    stop("'stat$gene' must match the graph vertex names.", call. = FALSE)
  }

  invisible(TRUE)
}


# Check scalar tuning parameters for causal GRN orientation.
.check_causalgrn_params <- function(
    alpha, conservative, max_order, max_dist, evidence
) {
  if (
    length(alpha) != 1L || !is.numeric(alpha) ||
      !is.finite(alpha) || alpha <= 0 || alpha >= 1
  ) {
    stop("'alpha' must be a single finite number between 0 and 1.", call. = FALSE)
  }
  if (
    !is.logical(conservative) || length(conservative) != 1L ||
      is.na(conservative)
  ) {
    stop("'conservative' must be TRUE or FALSE.", call. = FALSE)
  }
  if (
    length(max_order) != 1L || !is.numeric(max_order) ||
      !is.finite(max_order) || max_order != round(max_order) ||
      !(max_order %in% c(1, 2))
  ) {
    stop("'max_order' must be 1 or 2.", call. = FALSE)
  }
  if (
    length(max_dist) != 1L || !is.numeric(max_dist) ||
      is.na(max_dist) || max_dist < 1 ||
      (!is.infinite(max_dist) && max_dist != round(max_dist))
  ) {
    stop("'max_dist' must be a positive integer or Inf.", call. = FALSE)
  }
  if (
    length(evidence) != 1L || !is.numeric(evidence) ||
      !is.finite(evidence) || evidence < 1 || evidence != round(evidence)
  ) {
    stop("'evidence' must be a single positive integer.", call. = FALSE)
  }

  invisible(TRUE)
}
