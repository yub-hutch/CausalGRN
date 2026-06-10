#' Predict KO Effect Using the Mean Perturbation Effect Baseline
#'
#' This function calculates the average perturbation effect (delta from wild-type)
#' from a training set and uses it to predict the effect for a test set.
#'
#' @param Y A numeric matrix of expression data (cells x genes).
#' @param group Named character vector of cell labels: 'WT' for wild-type cells,
#'   perturbed gene for perturbed cell.
#' @param train_kos A character vector of the training knockout group names.
#' @param test_kos A character vector of the test knockout group names.
#' @param wt_name A character string specifying the wild-type group name (default is "WT").
#'
#' @return A numeric matrix of predicted delta values (rows are test KOs, columns are genes).
#' @export
predict_mean_perturbation_effect <- function(Y, group, train_kos, test_kos, wt_name = "WT") {
  .check_perturbation_effect_inputs(Y = Y, group = group)

  # 1. Perform aggressive data checks
  stopifnot(wt_name %in% group)
  stopifnot(all(train_kos %in% group))
  stopifnot(all(test_kos %in% group))

  # 2. Calculate the mean wild-type expression vector
  wt_expressions <- colMeans(Y[group == wt_name, , drop = FALSE])

  # 3. Calculate the "Mean Delta" from the training data
  train_delta_matrix <- vapply(train_kos, function(ko) {
    colMeans(Y[group == ko, , drop = FALSE]) - wt_expressions
  }, FUN.VALUE = numeric(ncol(Y)))

  # Average all the delta vectors together
  mean_delta <- rowMeans(train_delta_matrix)

  # 4. Create the final output matrix of predicted deltas
  # Each row is the identical mean delta vector.
  pred_matrix <- t(replicate(length(test_kos), mean_delta))
  rownames(pred_matrix) <- test_kos

  return(pred_matrix)
}
