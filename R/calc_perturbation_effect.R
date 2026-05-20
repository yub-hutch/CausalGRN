#' Calculate Perturbation Effect
#'
#' Calculate differential expression statistics between wild-type cells and
#' perturbed cells for each perturbation.
#'
#' @param Y Matrix of normalized scRNA-seq data of wild-type and perturbed cells.
#' @param group Named character vector of cell label: 'WT' for wild-type cells,
#' perturbed gene for perturbed cell.
#' @param ncores Number of CPUs to use.
#' @param min_cells Minimum number of cells required in each group.
#' @param gene_block_size Optional number of genes per block. Under the default
#' \code{NULL}, block size is chosen from available memory. Set this manually to
#' override automatic memory-aware blocking.
#' @return Tibble with columns:
#' \itemize{
#'   \item \code{ko}: Perturbed gene.
#'   \item \code{gene}: Affected gene.
#'   \item \code{diff}: Mean difference (Perturbed - WT).
#'   \item \code{cd}: Cohen's D using the pooled within-group SD.
#'   \item \code{wilcox_pv}: Wilcoxon rank sum test P-value.
#'   \item \code{t_pv}: T-test P-value.
#'   \item \code{cor_pearson}: Pearson correlation in wild-type cells.
#'   \item \code{cor_spearman}: Spearman correlation in wild-type cells.
#'   \item \code{wilcox_adj_pv}: BH-adjusted Wilcoxon rank sum test P-value.
#'   \item \code{t_adj_pv}: BH-adjusted T-test P-value.
#'   \item \code{adj_pv}: The larger of \code{wilcox_adj_pv} and \code{t_adj_pv}.
#' }
#'
#' @examples
#' # x -> y -> z
#' set.seed(123)
#' n <- 1e3
#' x0 <- rnorm(n)
#' x1 <- rnorm(n, mean = -1)
#' x <- c(x0, x1)
#' group <- c(rep('WT', n), rep('x', n))
#' y <- rnorm(2 * n, mean = x, sd = 1)
#' z <- rnorm(2 * n, mean = y, sd = 1)
#' data <- cbind(x = x, y = y, z = z)
#' rownames(data) <- paste0("cell", seq_len(nrow(data)))
#' names(group) <- rownames(data)
#' stat <- calc_perturbation_effect(
#'   Y = data,
#'   group = group,
#'   ncores = 2
#' )
#' print(stat)
#' @export
calc_perturbation_effect <- function(
    Y, group, ncores, min_cells = 50L, gene_block_size = NULL
) {
  .check_expression_matrix(Y)
  .check_group(group, min_cells = min_cells)
  .check_perturbation_effect_inputs(
    Y = Y,
    group = group
  )
  .check_ncores(ncores)
  .check_perturbation_effect_params(gene_block_size = gene_block_size)

  kos <- setdiff(group, 'WT')
  genes <- colnames(Y)
  wt_idx <- which(group == 'WT')
  ko_indices <- setNames(lapply(kos, \(ko) which(group == ko)), kos)

  njobs <- min(ncores, length(kos))
  ko_batch_size <- ceiling(length(kos) / njobs)
  ko_batches <- split(kos, ceiling(seq_along(kos) / ko_batch_size))

  if (is.null(gene_block_size)) {
    available_bytes <- .available_memory_bytes()
    if (is.na(available_bytes)) {
      block_size <- length(genes)
      message("Available memory could not be detected; using one gene batch.")
    } else {
      budget_fraction <- .memory_budget_fraction(available_bytes)
      budget_bytes <- available_bytes * budget_fraction
      rows_per_worker <- length(wt_idx) + max(vapply(
        ko_batches,
        \(batch) sum(lengths(ko_indices[batch])),
        numeric(1)
      ))
      bytes_per_gene <- 8 * rows_per_worker * 4
      block_size <- floor(budget_bytes / bytes_per_gene)
      block_size <- max(1L, min(length(genes), block_size))

      message(
        "Auto gene_block_size = ", block_size,
        " genes (available memory: ",
        round(available_bytes / 1024^3, 1),
        " GB; budget: ", budget_fraction * 100, "%; overhead: 4x)."
      )
    }
  } else {
    block_size <- gene_block_size
    message("Using gene_block_size = ", block_size, " genes.")
  }
  gene_blocks <- split(seq_along(genes), ceiling(seq_along(genes) / block_size))

  message(
    "Processing ", length(gene_blocks), " gene batch(es); ",
    length(ko_batches), " KO batch(es) per gene batch."
  )

  stat_by_ko <- setNames(vector('list', length(kos)), kos)
  for (ko in kos) {
    stat_by_ko[[ko]] <- vector('list', length(gene_blocks))
  }

  for (block_id in seq_along(gene_blocks)) {
    cols <- gene_blocks[[block_id]]
    wt_block <- Y[wt_idx, cols, drop = FALSE]
    block_genes <- genes[cols]
    message(
      "Gene batch ", block_id, "/", length(gene_blocks),
      ": ", length(cols), " genes."
    )

    # Workers receive only block-level slices, never the full Y matrix.
    jobs <- lapply(ko_batches, \(ko_batch) {
      pt_blocks <- setNames(
        lapply(ko_batch, \(ko) Y[ko_indices[[ko]], cols, drop = FALSE]),
        ko_batch
      )
      list(
        ko_batch = ko_batch,
        genes = block_genes,
        wt_block = wt_block,
        ko_wt_block = Y[wt_idx, ko_batch, drop = FALSE],
        pt_blocks = pt_blocks
      )
    })

    job_results <- .parallel_lapply(
      jobs,
      .calc_perturbation_effect_job,
      ncores = length(jobs),
      export = c(".calc_perturbation_effect_job", ".calc_cohens_d")
    )

    for (job_result in job_results) {
      for (ko in names(job_result)) {
        stat_by_ko[[ko]][[block_id]] <- job_result[[ko]]
      }
    }
  }

  stat <- do.call(
    rbind,
    unlist(stat_by_ko, recursive = FALSE, use.names = FALSE)
  )

  stat$wilcox_adj_pv <- p.adjust(stat$wilcox_pv, method = 'BH')
  stat$t_adj_pv <- p.adjust(stat$t_pv, method = 'BH')
  stat$adj_pv <- pmax(stat$wilcox_adj_pv, stat$t_adj_pv)
  return(stat)
}


# Compute perturbation-effect statistics for one KO batch and one gene block.
.calc_perturbation_effect_job <- function(job) {
  stat_list <- lapply(job$ko_batch, \(ko) {
    wt <- job$wt_block
    pt <- job$pt_blocks[[ko]]
    cd <- .calc_cohens_d(wt = wt, pt = pt)
    wilcox_pvs <- matrixTests::col_wilcoxon_twosample(
      wt,
      pt,
      exact = FALSE,
      correct = TRUE
    )$pvalue
    t_pvs <- matrixTests::col_t_welch(wt, pt)$pvalue
    cors_pearson <- c(cor(job$ko_wt_block[, ko], wt))
    cors_spearman <- c(cor(job$ko_wt_block[, ko], wt, method = 'spearman'))

    dplyr::tibble(
      ko = ko, gene = job$genes,
      diff = cd$diff, cd = cd$cd,
      wilcox_pv = wilcox_pvs, t_pv = t_pvs,
      cor_pearson = cors_pearson, cor_spearman = cors_spearman
    )
  })
  names(stat_list) <- job$ko_batch
  return(stat_list)
}


# Read available physical memory in bytes. Return NA when unavailable.
.available_memory_bytes <- function() {
  sysname <- Sys.info()[["sysname"]]

  if (identical(sysname, "Linux")) {
    if (!file.exists("/proc/meminfo")) {
      return(NA_real_)
    }

    meminfo <- readLines("/proc/meminfo", warn = FALSE)
    line <- grep("^MemAvailable:", meminfo, value = TRUE)
    if (!length(line)) {
      return(NA_real_)
    }

    kb <- suppressWarnings(as.numeric(sub(
      "^MemAvailable:\\s+([0-9]+).*",
      "\\1",
      line[1]
    )))
    if (is.na(kb) || kb <= 0) {
      return(NA_real_)
    }

    return(kb * 1024)
  }

  if (identical(sysname, "Darwin")) {
    page_size <- suppressWarnings(as.numeric(system2(
      "sysctl",
      c("-n", "hw.pagesize"),
      stdout = TRUE,
      stderr = FALSE
    )[1]))
    vm <- suppressWarnings(system2(
      "vm_stat",
      stdout = TRUE,
      stderr = FALSE
    ))
    if (is.na(page_size) || page_size <= 0 || !length(vm)) {
      return(NA_real_)
    }

    get_pages <- function(label) {
      line <- grep(paste0("^", label, ":"), vm, value = TRUE)
      if (!length(line)) {
        return(0)
      }
      value <- gsub("[^0-9]", "", line[1])
      suppressWarnings(as.numeric(value))
    }

    pages <- sum(
      get_pages("Pages free"),
      get_pages("Pages inactive"),
      get_pages("Pages speculative"),
      na.rm = TRUE
    )
    if (is.na(pages) || pages <= 0) {
      return(NA_real_)
    }

    return(pages * page_size)
  }

  if (identical(.Platform$OS.type, "windows")) {
    kb <- suppressWarnings(as.numeric(system2(
      "powershell",
      c(
        "-NoProfile",
        "-Command",
        "(Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory"
      ),
      stdout = TRUE,
      stderr = FALSE
    )[1]))
    if (is.na(kb) || kb <= 0) {
      wmic <- suppressWarnings(system2(
        "wmic",
        c("OS", "get", "FreePhysicalMemory", "/Value"),
        stdout = TRUE,
        stderr = FALSE
      ))
      line <- grep("^FreePhysicalMemory=", wmic, value = TRUE)
      if (length(line)) {
        kb <- suppressWarnings(as.numeric(sub(
          "^FreePhysicalMemory=",
          "",
          line[1]
        )))
      }
    }
    if (is.na(kb) || kb <= 0) {
      return(NA_real_)
    }

    return(kb * 1024)
  }

  NA_real_
}


# Use more of available memory on large-memory machines.
.memory_budget_fraction <- function(available_bytes) {
  gb <- available_bytes / 1024^3
  if (gb < 20) {
    return(0.25)
  }
  if (gb < 100) {
    return(0.50)
  }
  if (gb < 500) {
    return(0.70)
  }

  0.80
}


# Standard two-sample Cohen's d for perturbed versus wild-type cells.
.calc_cohens_d <- function(wt, pt) {
  diffs <- colMeans(pt) - colMeans(wt)
  wt_sds <- apply(wt, 2, sd)
  pt_sds <- apply(pt, 2, sd)
  pooled_sds <- sqrt(
    ((nrow(wt) - 1) * wt_sds ^ 2 + (nrow(pt) - 1) * pt_sds ^ 2) /
      (nrow(wt) + nrow(pt) - 2)
  )

  return(list(diff = diffs, cd = diffs / pooled_sds))
}
