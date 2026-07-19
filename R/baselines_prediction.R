#' Predict KO Effect Using the Mean Perturbation Effect Baseline
#'
#' Calculates the average perturbation effect (delta from wild-type) from a
#' training set and uses it to predict the effect for a test set.
#'
#' @param Y A numeric matrix of expression data (cells x genes).
#' @param group Named character vector of cell labels: 'WT' for wild-type cells,
#'   perturbed gene for perturbed cell.
#' @param train_kos A character vector of the training knockout group names.
#' @param test_kos A character vector of the test knockout group names.
#'
#' @return A numeric matrix of predicted delta values (rows are test KOs,
#'   columns are genes).
#' @export
predict_mean_perturbation_effect <- function(
    Y, group, train_kos, test_kos
) {
  .check_perturbation_effect_inputs(Y = Y, group = group)

  if (
    !is.character(train_kos) || length(train_kos) == 0L ||
      anyNA(train_kos) || any(train_kos == '')
  ) {
    stop("'train_kos' must be a non-empty character vector.", call. = FALSE)
  }
  if (
    !is.character(test_kos) || length(test_kos) == 0L ||
      anyNA(test_kos) || any(test_kos == '')
  ) {
    stop("'test_kos' must be a non-empty character vector.", call. = FALSE)
  }
  if (!all(train_kos %in% group)) {
    stop("'train_kos' must all be present in 'group'.", call. = FALSE)
  }
  if (!all(test_kos %in% group)) {
    stop("'test_kos' must all be present in 'group'.", call. = FALSE)
  }

  wt_expressions <- colMeans(Y[group == 'WT', , drop = FALSE])

  train_delta_matrix <- vapply(
    train_kos,
    function(ko) {
      colMeans(Y[group == ko, , drop = FALSE]) - wt_expressions
    },
    FUN.VALUE = numeric(ncol(Y))
  )

  mean_delta <- rowMeans(train_delta_matrix)

  pred_matrix <- matrix(
    mean_delta,
    nrow = length(test_kos),
    ncol = length(mean_delta),
    byrow = TRUE
  )
  colnames(pred_matrix) <- colnames(Y)
  rownames(pred_matrix) <- test_kos
  return(pred_matrix)
}


#' Predict Oracle Perturbation Effect
#'
#' Predicts the perturbation effect (delta from wild-type) for test cells using
#' a fitted B matrix, where the true expression values of the predictors are
#' known.
#'
#' @param B A numeric matrix of regulatory coefficients where rows are sources
#'   ('Intercept' and genes) and columns are targets (genes). This is the output
#'   of \code{\link{fit_expression_model}}.
#' @param Y A numeric matrix of expression data (cells x genes) for which to make
#'   predictions.
#' @param group Named character vector indicating the group for each cell in
#'   \code{Y}. Wild-type cells are not required for oracle prediction.
#' @param wt_expressions A named numeric vector of the wild-type expression
#'   levels for each gene.
#' @return A numeric matrix of predicted delta values (rows are unique KOs,
#'   columns are genes).
#' @export
predict_oracle_perturbation_effect <- function(
    B, Y, group, wt_expressions
) {
  .check_expression_matrix(Y)

  genes <- colnames(Y)
  .check_expression_model_matrix(B, genes = genes)
  .check_expression_vector(
    wt_expressions,
    arg = 'wt_expressions',
    expected_names = genes
  )
  .check_group(group, require_wt = FALSE)
  if (length(group) != nrow(Y)) {
    stop("'group' must have one label per row of 'Y'.", call. = FALSE)
  }
  if (!identical(names(group), rownames(Y))) {
    stop("'group' names must be identical to 'Y' row names.", call. = FALSE)
  }

  X_full <- cbind(Intercept = 1, Y)
  pred_cell_abs <- X_full %*% B

  kos <- sort(unique(group))
  pred_ko_abs <- matrix(
    vapply(
      kos,
      function(ko) {
        colMeans(pred_cell_abs[group == ko, , drop = FALSE])
      },
      FUN.VALUE = numeric(length(genes))
    ),
    nrow = length(kos),
    ncol = length(genes),
    byrow = TRUE
  )
  colnames(pred_ko_abs) <- genes
  rownames(pred_ko_abs) <- kos

  pred_delta <- pred_ko_abs - matrix(
    wt_expressions,
    nrow = nrow(pred_ko_abs),
    ncol = ncol(pred_ko_abs),
    byrow = TRUE
  )
  return(pred_delta)
}
