sample_BA_dag <- function(d, nodes) {
  raw <- pcalg::randDAG(n = length(nodes), d = d, method = 'barabasi', DAG = TRUE, weighted = FALSE)
  graph <- igraph::graph_from_adjacency_matrix(as(raw, 'matrix'), mode = 'directed')
  igraph::V(graph)$name <- nodes
  return(graph)
}


simulate_coef <- function(parents, b0, min_coef, max_coef) {
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
#' @param center_normal_sd Center of standard deviations of latent normal variables. Gene specific SDs are sampled from Gamma distribution.
#' @param center_ko_eff Center of KO effectiveness. KO specific effectiveness is sampled from Beta distribution.
#' @return List containing:
#' \itemize{
#'   \item \code{dag}: Simulated DAG.
#'   \item \code{coef}: Simulated coefficients.
#'   \item \code{Y}: Simulated expression.
#'   \item \code{group}: Simulated group.
#' }
#' @export
simulate_grn_guided_expression <- function(d, Y, group, min_coef = 0.3, max_coef = 0.5, center_normal_sd = 1, center_ko_eff = 0.9) {
  # Check inputs
  stopifnot(identical(rownames(Y), names(group)))
  stopifnot('WT' %in% group)
  stopifnot(all(table(group) >= 50))
  stopifnot(center_normal_sd > 0)
  stopifnot(center_ko_eff > 0 & center_ko_eff < 1)
  kos <- setdiff(group, 'WT')
  genes <- colnames(Y)
  stopifnot(all(kos %in% genes))
  ngene <- length(genes)
  # Simulate DAG
  dag <- sample_BA_dag(d = d, nodes = genes)
  ordering <- igraph::topo_sort(dag, mode = 'out')$name
  # Initialize expression matrix
  nwt <- sum(group == 'WT')
  wt <- u_wt <- matrix(
    data = NA,
    nrow = nwt,
    ncol = ngene,
    dimnames = list(names(group)[group == 'WT'], genes)
  )
  npt <- sapply(kos, \(ko) sum(group == ko))
  pts <- u_pts <- sapply(kos, simplify = F, \(ko) {
    matrix(
      data = NA,
      nrow = npt[ko],
      ncol = ngene,
      dimnames = list(names(group)[group == ko], genes)
    )
  })
  # Set parameters
  max_wt_value <- max(Y[group == 'WT', ])
  wt_colmeans <- colMeans(Y[group == 'WT', ])
  max_wt_colmean <- max(wt_colmeans)
  lib_sizes <- rowSums(Y)
  lib_sizes <- lib_sizes / median(lib_sizes)
  normal_sds <- sqrt(rgamma(n = ngene, shape = 2, rate = 2 / (center_normal_sd ^ 2)))
  names(normal_sds) <- genes
  # Simulate expression
  coef <- setNames(vector('list', ngene), genes)
  for (v in ordering) {
    parents <- igraph::neighbors(graph = dag, v = v, mode = 'in')$name
    if (length(parents) == 0) {
      u_wt[, v] <- rnorm(n = nwt, mean = log(wt_colmeans[v]), sd = normal_sds[v])
      wt[, v] <- rpois(n = nwt, lambda = lib_sizes[group == 'WT'] * exp(u_wt[, v]))
      wt[, v] <- pmin(wt[, v], max_wt_value)
      for (ko in kos) {
        u_pts[[ko]][, v] <- rnorm(n = npt[ko], mean = log(wt_colmeans[v]), sd = normal_sds[v])
        if (ko == v) {
          eff <- rbeta(n = npt[ko], shape1 = center_ko_eff * 50, shape2 = (1 - center_ko_eff) * 50)
          u_pts[[ko]][, v] <- u_pts[[ko]][, v] + log(1 - eff)
        }
        pts[[ko]][, v] <- rpois(n = npt[ko], lambda = lib_sizes[group == ko] * exp(u_pts[[ko]][, v]))
        pts[[ko]][, v] <- pmin(pts[[ko]][, v], max_wt_value)
      }
    } else {
      while (TRUE) {
        coef[[v]] <- simulate_coef(parents, b0 = log(wt_colmeans[v]), min_coef = min_coef, max_coef = max_coef)
        normal_means <- as.vector(coef[[v]][1] + u_wt[, parents, drop = FALSE] %*% coef[[v]][-1])
        u_wt[, v] <- rnorm(n = nwt, mean = normal_means, sd = normal_sds[v])
        wt[, v] <- rpois(n = nwt, lambda = lib_sizes[group == 'WT'] * exp(u_wt[, v]))
        wt[, v] <- pmin(wt[, v], max_wt_value)
        if (mean(wt[, v]) <= max_wt_colmean) break # Make sure the sampled coefficients are reasonable
      }
      for (ko in kos) {
        normal_means <- as.vector(coef[[v]][1] + u_pts[[ko]][, parents, drop = FALSE] %*% coef[[v]][-1])
        u_pts[[ko]][, v] <- rnorm(n = npt[ko], mean = normal_means, sd = normal_sds[v])
        if (ko == v) {
          eff <- rbeta(n = npt[ko], shape1 = center_ko_eff * 50, shape2 = (1 - center_ko_eff) * 50)
          u_pts[[ko]][, v] <- u_pts[[ko]][, v] + log(1 - eff)
        }
        pts[[ko]][, v] <- rpois(n = npt[ko], lambda = lib_sizes[group == ko] * exp(u_pts[[ko]][, v]))
        pts[[ko]][, v] <- pmin(pts[[ko]][, v], max_wt_value)
      }
    }
  }
  Y <- do.call(rbind, c(list(wt), pts))
  group <- group[rownames(Y)]
  return(list(dag = dag, coef = coef, Y = Y, group = group))
}
