sample_BA_dag <- function(d, nodes) {
  raw <- pcalg::randDAG(n = length(nodes), d = d, method = 'barabasi', DAG = TRUE, weighted = FALSE)
  graph <- igraph::graph_from_adjacency_matrix(as(raw, 'matrix'), mode = 'directed')
  igraph::V(graph)$name <- nodes
  return(graph)
}


simulate_coef <- function(parents, real_colmeans, min_coef, max_coef) {
  stopifnot(length(parents) > 0)
  b0 <- sample(x = log(real_colmeans), size = 1)
  b1 <- runif(n = length(parents), min = min_coef, max = max_coef)
  b1 <- b1 * sample(x = c(-1, 1), size = length(b1), replace = TRUE)
  names(b1) <- parents
  return(c(Intercept = as.vector(b0), b1))
}


#' GRN-guided simulation of scRNA-seq data of wild-type and perturbed cells mimicking real data
#'
#' Simulates gene expression based on simulated directed acyclic graph and real single-cell perturbation data.
#'
#' @param d Expected number of neighbors (in/out).
#' @param Y Real single-cell perturbation data (cells x genes).
#' @param group Named vector of cell label: 'WT' for wild-type cells, perturbed gene for perturbed cell.
#' @param min_coef Minimum value of simulated coefficients.
#' @param max_coef Maximum value of simulated coefficients.
#' @param cv Coefficient of variation for simulated normal distributions.
#' @param lib_sizes Vector of cell library sizes.
#' @return List containing:
#' \itemize{
#'   \item \code{dag}: Simulated DAG.
#'   \item \code{coef}: Simulated coefficients.
#'   \item \code{Y}: Simulated expression.
#'   \item \code{group}: Simulated group.
#' }
#' @export
simulate_grn_guided_expression <- function(d, Y, group, min_coef, max_coef, cv, lib_sizes) {
  # Check inputs
  stopifnot(identical(rownames(Y), names(group)))
  stopifnot('WT' %in% group)
  stopifnot(all(table(group) >= 50))
  stopifnot(length(lib_sizes) == nrow(Y))
  kos <- setdiff(group, 'WT')
  genes <- colnames(Y)
  stopifnot(all(kos %in% genes))
  nko <- length(kos)
  ngene <- length(genes)
  # Simulate DAG
  dag <- sample_BA_dag(d = d, nodes = genes)
  ordering <- igraph::topo_sort(dag, mode = 'out')$name
  # Initialize expression matrix
  nwt <- sum(group == 'WT')
  wt <- matrix(
    data = NA,
    nrow = nwt,
    ncol = ngene,
    dimnames = list(names(group)[group == 'WT'], genes)
  )
  npt <- sapply(kos, \(ko) sum(group == ko))
  pts <- sapply(kos, simplify = F, \(ko) {
    matrix(
      data = NA,
      nrow = npt[ko],
      ncol = ngene,
      dimnames = list(names(group)[group == ko], genes)
    )
  })
  coef <- setNames(vector('list', ngene), genes)
  # Set empirical rules
  real_wt <- Y[group == 'WT', ]
  max_wt_value <- max(real_wt)
  wt_colmeans <- colMeans(real_wt)
  max_wt_colmean <- max(wt_colmeans)
  # Simulate expression
  pb <- pbmcapply::progressBar(min = 0, max = ngene)
  for (v in ordering) {
    setTxtProgressBar(pb = pb, value = match(v, ordering))
    parents <- igraph::neighbors(graph = dag, v = v, mode = 'in')$name
    if (length(parents) == 0) {
      wt[, v] <- as.integer(sample(x = real_wt[, v]) * lib_sizes[group == 'WT'])
      for (ko in kos) {
        if (ko == v) {
          pts[[ko]][, v] <- as.integer(sample(x = Y[group == ko, v]) * lib_sizes[group == ko])
        } else {
          pts[[ko]][, v] <- as.integer(sample(x = real_wt[, v], size = npt[ko], replace = TRUE) * lib_sizes[group == ko])
        }
      }
    } else {
      accept <- FALSE
      while (!accept) {
        coef[[v]] <- simulate_coef(
          parents = parents,
          real_colmeans = wt_colmeans,
          min_coef = min_coef,
          max_coef = max_coef
        )
        normal_means <- as.vector(coef[[v]][1] + log(1 + wt[, parents, drop = FALSE]) %*% coef[[v]][-1])
        normal_sds <- abs(normal_means / cv)
        lambdas <- exp(rnorm(n = nwt, mean = normal_means, sd = normal_sds))
        wt[, v] <- rpois(n = nwt, lambda = lambdas * lib_sizes[group == 'WT'])
        wt[, v] <- pmin(wt[, v], max_wt_value)
        accept <- (mean(wt[, v]) <= max_wt_colmean)
      }
      for (ko in kos) {
        if (ko == v) {
          pts[[ko]][, v] <- 0
        } else {
          normal_means <- as.vector(coef[[v]][1] + log(1 + pts[[ko]][, parents, drop = FALSE]) %*% coef[[v]][-1])
          normal_sds <- abs(normal_means / cv)
          lambdas <- exp(rnorm(n = npt[ko], mean = normal_means, sd = normal_sds))
          pts[[ko]][, v] <- rpois(n = npt[ko], lambda = lambdas * lib_sizes[group == ko])
          pts[[ko]][, v] <- pmin(pts[[ko]][, v], max_wt_value)
        }
      }
    }
  }
  close(pb)
  Y <- do.call(rbind, c(list(wt), pts))
  group <- group[rownames(Y)]
  return(list(dag = dag, coef = coef, Y = Y, group = group))
}
