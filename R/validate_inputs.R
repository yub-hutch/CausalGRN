# Check group labels and return them as a named character vector.
.check_group <- function(group, row_names = NULL, require_wt = TRUE, min_cells = NULL) {
  if (!is.atomic(group) || !is.null(dim(group))) {
    stop("'group' must be an atomic vector.", call. = FALSE)
  }
  if (is.null(names(group))) {
    stop("'group' must be named.", call. = FALSE)
  }

  group_names <- names(group)
  group <- as.character(group)
  names(group) <- group_names

  if (any(is.na(group) | group == "")) {
    stop("'group' must not contain missing or empty labels.", call. = FALSE)
  }
  if (!is.null(row_names) && !identical(row_names, names(group))) {
    stop("'group' names must be identical to matrix row names.", call. = FALSE)
  }
  if (require_wt && !"WT" %in% group) {
    stop("'group' must contain wild-type cells labeled as 'WT'.", call. = FALSE)
  }
  if (!is.null(min_cells) && any(table(group) < min_cells)) {
    stop("Each group must contain at least ", min_cells, " cells.", call. = FALSE)
  }

  return(group)
}


# Check core count input and return ncores as a positive integer.
.check_ncores <- function(ncores) {
  if (length(ncores) != 1L || is.na(ncores) || !is.finite(ncores)) {
    stop("'ncores' must be a single finite number.", call. = FALSE)
  }

  ncores <- as.integer(ncores)
  if (ncores < 1L) {
    stop("'ncores' must be >= 1.", call. = FALSE)
  }

  return(ncores)
}


# Check that count data are a dense, finite, non-negative, integer-like matrix.
.check_count_matrix <- function(count, group = NULL) {
  if (!is.matrix(count) || !is.numeric(count)) {
    stop("'count' must be a numeric matrix.", call. = FALSE)
  }
  if (is.null(rownames(count)) || is.null(colnames(count))) {
    stop("'count' must have row and column names.", call. = FALSE)
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

  if (!is.null(group)) {
    .check_group(group, row_names = rownames(count), require_wt = FALSE)
  }

  return(count)
}


# Check that normalized expression data are a finite numeric matrix.
.check_expression_matrix <- function(Y, group = NULL) {
  if (!is.matrix(Y) || !is.numeric(Y)) {
    stop("'Y' must be a numeric matrix.", call. = FALSE)
  }
  if (is.null(rownames(Y)) || is.null(colnames(Y))) {
    stop("'Y' must have row and column names.", call. = FALSE)
  }

  for (j in seq_len(ncol(Y))) {
    if (any(!is.finite(Y[, j]))) {
      stop("'Y' must contain only finite values.", call. = FALSE)
    }
  }

  if (!is.null(group)) {
    .check_group(group, row_names = rownames(Y), require_wt = FALSE)
  }

  return(Y)
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
