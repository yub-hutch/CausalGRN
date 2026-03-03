.causalgrn_normalize_ncores <- function(ncores) {
  if (length(ncores) != 1L || is.na(ncores) || !is.finite(ncores)) {
    stop("'ncores' must be a single, finite number.")
  }
  ncores <- as.integer(ncores)
  if (ncores < 1L) {
    stop("'ncores' must be >= 1.")
  }
  ncores
}

.causalgrn_parallel_lapply <- function(
    X,
    FUN,
    ...,
    ncores = 1,
    preschedule = TRUE,
    export = NULL
) {
  ncores <- .causalgrn_normalize_ncores(ncores)

  if (ncores <= 1L) {
    return(lapply(X, FUN, ...))
  }

  if (.Platform$OS.type != "windows") {
    return(pbmcapply::pbmclapply(
      X,
      FUN,
      ...,
      mc.cores = ncores,
      mc.preschedule = preschedule
    ))
  }

  cl <- parallel::makeCluster(ncores)
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)

  parallel::clusterEvalQ(cl, {
    if (!requireNamespace("CausalGRN", quietly = TRUE)) {
      stop("Package 'CausalGRN' is required on PSOCK workers.")
    }
    NULL
  })

  if (!is.null(export) && length(export) > 0L) {
    parallel::clusterExport(cl, varlist = export, envir = parent.frame())
  }

  old_pboptions <- pbapply::pboptions(type = "timer", use_lb = isFALSE(preschedule))
  on.exit(pbapply::pboptions(old_pboptions), add = TRUE)

  tryCatch(
    pbapply::pblapply(X, FUN, ..., cl = cl),
    interrupt = function(e) {
      try(parallel::stopCluster(cl), silent = TRUE)
      stop(e)
    }
  )
}
