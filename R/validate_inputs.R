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

  group
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

  count
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

  Y
}
