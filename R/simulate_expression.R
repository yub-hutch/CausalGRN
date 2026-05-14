# Sample a Barabasi-Albert DAG over the provided nodes.
.sample_ba_dag <- function(d, nodes) {
  raw <- pcalg::randDAG(
    n = length(nodes),
    d = d,
    method = 'barabasi',
    DAG = TRUE,
    weighted = FALSE
  )
  graph <- igraph::graph_from_adjacency_matrix(as(raw, 'matrix'), mode = 'directed')
  igraph::V(graph)$name <- nodes
  return(graph)
}


# Sample signed parent coefficients and combine them with the intercept.
.simulate_coef <- function(parents, b0, min_coef, max_coef) {
  b1 <- runif(n = length(parents), min = min_coef, max = max_coef)
  b1 <- b1 * sample(x = c(-1, 1), size = length(b1), replace = TRUE)
  names(b1) <- parents
  return(c(Intercept = as.vector(b0), b1))
}

#' Simulate GRN-guided count data mimicking real single-cell perturbation data
#'
#' Simulates count data for wild-type and perturbed cells using a synthetic
#' directed acyclic graph and a reference count matrix to mimic real data.
#'
#' @param d Expected total degree per gene (parents plus children) in the DAG.
#' @param count Reference count matrix (cells x genes).
#' @param group Named vector of cell labels. Use \code{"WT"} for wild-type cells
#'   and gene names for perturbed cells.
#' @param min_coef Minimum absolute regulatory coefficient.
#' @param max_coef Maximum absolute regulatory coefficient.
#' @param center_normal_sd Center of latent normal standard deviations. Gene-specific
#'   variances are sampled from a Gamma distribution.
#' @param center_ko_eff Center of knockout efficacy. Cell-specific efficacies are
#'   sampled from a Beta distribution.
#' @param max_attempts Maximum attempts to sample stable coefficients per gene.
#' @return List containing:
#' \itemize{
#'   \item \code{dag}: Simulated DAG.
#'   \item \code{coef}: Simulated coefficients.
#'   \item \code{count}: Simulated count matrix.
#'   \item \code{group}: Simulated group.
#' }
#' @examples
#' # Create a small reference count matrix mimicking real perturbation data.
#' set.seed(1)
#' genes <- paste0("Gene", 1:5)
#' group <- rep(c("WT", "Gene1"), each = 50)
#' names(group) <- paste0("cell", seq_along(group))
#' count <- matrix(
#'   rpois(length(group) * length(genes), lambda = 10),
#'   nrow = length(group),
#'   dimnames = list(names(group), genes)
#' )
#'
#' # Simulate count data from a synthetic GRN.
#' sim <- simulate_grn_guided_expression(d = 2, count = count, group = group)
#' str(sim)
#' @export
simulate_grn_guided_expression <- function(
    d, count, group, min_coef = 0.3, max_coef = 0.5,
    center_normal_sd = 1, center_ko_eff = 0.9, max_attempts = 10000
) {
  # Validate inputs
  count <- .check_count_matrix(count)
  group <- .check_group(group, row_names = rownames(count), require_wt = TRUE, min_cells = 50)

  if (length(d) != 1 || !is.finite(d) || d <= 1) {
    stop("'d' must be a single finite number > 1.", call. = FALSE)
  }

  if (length(min_coef) != 1 || !is.finite(min_coef) || min_coef < 0) {
    stop("'min_coef' must be a single finite non-negative number.", call. = FALSE)
  }

  if (length(max_coef) != 1 || !is.finite(max_coef) || max_coef < min_coef) {
    stop("'max_coef' must be a single finite number >= 'min_coef'.", call. = FALSE)
  }

  if (length(center_normal_sd) != 1 || !is.finite(center_normal_sd) || center_normal_sd <= 0) {
    stop("'center_normal_sd' must be a single finite positive number.", call. = FALSE)
  }

  if (length(center_ko_eff) != 1 || !is.finite(center_ko_eff) || center_ko_eff <= 0 || center_ko_eff >= 1) {
    stop("'center_ko_eff' must be a single finite number between 0 and 1.", call. = FALSE)
  }

  if (length(max_attempts) != 1 || !is.finite(max_attempts) || max_attempts < 1 || max_attempts != round(max_attempts)) {
    stop("'max_attempts' must be a single positive integer.", call. = FALSE)
  }

  kos <- setdiff(group, 'WT')
  genes <- colnames(count)
  invalid_kos <- setdiff(kos, genes)
  if (length(invalid_kos) > 0) {
    stop(
      "Perturbation labels must be gene names in 'count': ",
      paste(utils::head(invalid_kos, 10), collapse = ", "),
      if (length(invalid_kos) > 10) ", ..." else "",
      call. = FALSE
    )
  }

  ngene <- length(genes)
  if (d > ngene - 1) {
    stop("'d' must be <= the number of genes - 1.", call. = FALSE)
  }

  # Simulate DAG
  dag <- .sample_ba_dag(d = d, nodes = genes)
  ordering <- igraph::topo_sort(dag, mode = 'out')$name

  # Initialize simulated matrices
  nwt <- sum(group == 'WT')
  wt <- u_wt <- matrix(
    data = NA,
    nrow = nwt,
    ncol = ngene,
    dimnames = list(names(group)[group == 'WT'], genes)
  )
  npt <- vapply(kos, \(ko) sum(group == ko), numeric(1))
  pts <- u_pts <- setNames(lapply(kos, \(ko) {
    matrix(
      data = NA,
      nrow = npt[ko],
      ncol = ngene,
      dimnames = list(names(group)[group == ko], genes)
    )
  }), kos)

  # Set parameters
  max_wt_value <- max(count[group == 'WT', ])
  wt_colmeans <- colMeans(count[group == 'WT', ])
  zero_wt_mean_genes <- genes[wt_colmeans == 0]
  if (length(zero_wt_mean_genes) > 0) {
    stop(
      "Remove genes with zero WT mean before simulation: ",
      paste(utils::head(zero_wt_mean_genes, 10), collapse = ", "),
      if (length(zero_wt_mean_genes) > 10) ", ..." else "",
      call. = FALSE
    )
  }

  max_wt_colmean <- max(wt_colmeans)
  lib_sizes <- rowSums(count)
  lib_sizes <- lib_sizes / median(lib_sizes)
  normal_sds <- sqrt(rgamma(n = ngene, shape = 2, rate = 2 / (center_normal_sd ^ 2)))
  names(normal_sds) <- genes

  # Simulate count data
  coef <- setNames(vector('list', ngene), genes)
  for (v in ordering) {
    parents <- igraph::neighbors(graph = dag, v = v, mode = 'in')$name
    b0 <- log(wt_colmeans[v]) - 0.5 * normal_sds[v] ^ 2
    if (length(parents) == 0) {
      u_wt[, v] <- rnorm(n = nwt, mean = b0, sd = normal_sds[v])
      wt[, v] <- rpois(n = nwt, lambda = lib_sizes[group == 'WT'] * exp(u_wt[, v]))
      wt[, v] <- pmin(wt[, v], max_wt_value)
      for (ko in kos) {
        u_pts[[ko]][, v] <- rnorm(n = npt[ko], mean = b0, sd = normal_sds[v])
        if (ko == v) {
          eff <- rbeta(n = npt[ko], shape1 = center_ko_eff * 50, shape2 = (1 - center_ko_eff) * 50)
          u_pts[[ko]][, v] <- u_pts[[ko]][, v] + log(1 - eff)
        }
        pts[[ko]][, v] <- rpois(n = npt[ko], lambda = lib_sizes[group == ko] * exp(u_pts[[ko]][, v]))
        pts[[ko]][, v] <- pmin(pts[[ko]][, v], max_wt_value)
      }
    } else {
      for (attempt in seq_len(max_attempts)) {
        coef[[v]] <- .simulate_coef(parents, b0 = b0, min_coef = min_coef, max_coef = max_coef)
        normal_means <- as.vector(coef[[v]][1] + u_wt[, parents, drop = FALSE] %*% coef[[v]][-1])
        u_wt[, v] <- rnorm(n = nwt, mean = normal_means, sd = normal_sds[v])
        wt[, v] <- rpois(n = nwt, lambda = lib_sizes[group == 'WT'] * exp(u_wt[, v]))
        wt[, v] <- pmin(wt[, v], max_wt_value)
        if (mean(wt[, v]) <= max_wt_colmean) {
          break
        }
      }

      if (attempt == max_attempts && mean(wt[, v]) > max_wt_colmean) {
        stop(
          "Failed to sample stable coefficients for gene '", v, "' after ",
          max_attempts, " attempts. Increase max_attempts, reduce d or max_coef, ",
          "or filter genes with extreme WT expression.",
          call. = FALSE
        )
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

  count <- do.call(rbind, c(list(wt), pts))
  group <- group[rownames(count)]
  return(list(dag = dag, coef = coef, count = count, group = group))
}
