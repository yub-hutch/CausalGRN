#' Predict Oracle Perturbation Effect (Delta)
#'
#' Predicts the perturbation effect (delta from wild-type) for test cells using
#' a fitted B matrix, where the true expression values of the predictors are
#' known (oracle setting).
#'
#' @param B A numeric matrix of regulatory coefficients where rows are sources
#'   ('Intercept' and genes) and columns are targets (genes). This is the output
#'   of \code{\link{fit_expression_model}}.
#' @param test_Y A numeric matrix of expression data (cells x genes) for which to
#'   make predictions.
#' @param test_group A character or factor vector indicating the group for each
#'   cell in `test_Y`.
#' @param wt_expressions A named numeric vector of the wild-type expression levels
#'   for each gene.
#' @return A numeric matrix of predicted delta values (rows are unique KOs,
#'   columns are genes).
#' @export
predict_oracle_perturbation_effect <- function(B, test_Y, test_group, wt_expressions) {
  # --- 1. Input Validation ---
  .check_expression_matrix(test_Y, arg = "test_Y")
  if (!(is.character(test_group) || is.factor(test_group)) || !is.null(dim(test_group))) {
    stop("'test_group' must be a character or factor vector.", call. = FALSE)
  }
  if (length(test_group) != nrow(test_Y)) {
    stop("'test_group' must have one label per row of 'test_Y'.", call. = FALSE)
  }
  test_group <- as.character(test_group)
  if (any(is.na(test_group) | test_group == "")) {
    stop("'test_group' must not contain missing or empty labels.", call. = FALSE)
  }

  genes <- colnames(test_Y)
  stopifnot(
    is.matrix(B), is.matrix(test_Y),
    is.numeric(wt_expressions), !is.null(names(wt_expressions)),
    identical(colnames(B), genes),
    identical(rownames(B), c('Intercept', genes)),
    identical(names(wt_expressions), genes)
  )

  # --- 2. Perform Oracle Prediction for Absolute Expression ---
  # Prepare the full input matrix by adding an intercept column to the test data.
  X_full <- cbind(Intercept = 1, test_Y)

  # The core prediction is a single, efficient matrix multiplication.
  pred_cell_abs <- X_full %*% B

  # --- 3. Aggregate Cell-Level Predictions to KO-Level ---
  kos <- sort(unique(as.character(test_group)))

  # Aggregate by taking the mean absolute prediction for all cells in a KO group.
  pred_ko_abs <- vapply(kos, function(ko) {
    colMeans(pred_cell_abs[test_group == ko, , drop = FALSE])
  }, FUN.VALUE = numeric(length(genes)))

  pred_ko_abs <- t(pred_ko_abs)

  # --- 4. Convert Absolute Predictions to Delta Predictions ---
  # Subtract the WT expression from each row to get the predicted effect.
  pred_delta <- pred_ko_abs - matrix(wt_expressions, nrow = nrow(pred_ko_abs), ncol = ncol(pred_ko_abs), byrow = TRUE)

  return(pred_delta)
}
