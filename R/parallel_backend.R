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
  FUN <- match.fun(FUN)
  nout <- max(1L, min(1000L, length(X)))

  if (ncores <= 1L) {
    if (!interactive()) {
      return(lapply(X, FUN, ...))
    }

    old_pboptions <- pbapply::pboptions(type = "timer", nout = nout)
    on.exit(pbapply::pboptions(old_pboptions), add = TRUE)

    return(pbapply::pblapply(X, FUN, ...))
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
  pids <- try(unlist(parallel::clusterCall(cl, base::Sys.getpid)), silent = TRUE)

  .causalgrn_stop_cluster <- function(cl, pids = NULL, hard = FALSE) {
    for (node in cl) {
      try(parallel::sendData(node, list(type = "DONE", data = NULL, tag = NULL)), silent = TRUE)
      try(close(node$con), silent = TRUE)
    }
    if (isTRUE(hard) && is.numeric(pids)) {
      for (pid in pids) {
        try(tools::pskill(pid), silent = TRUE)
      }
    }
    invisible(TRUE)
  }

  on.exit(.causalgrn_stop_cluster(cl), add = TRUE)

  if (!is.primitive(FUN)) {
    environment(FUN) <- globalenv()
  }

  parallel::clusterEvalQ(cl, {
    if (!requireNamespace("CausalGRN", quietly = TRUE)) {
      stop("Package 'CausalGRN' is required on PSOCK workers.")
    }
    NULL
  })

  if (!is.null(export) && length(export) > 0L) {
    parallel::clusterExport(cl, varlist = export, envir = parent.frame())
  }

  old_pboptions <- pbapply::pboptions(
    type = "timer",
    use_lb = isFALSE(preschedule),
    nout = nout
  )
  on.exit(pbapply::pboptions(old_pboptions), add = TRUE)

  tryCatch(
    pbapply::pblapply(X, FUN, ..., cl = cl),
    interrupt = function(e) {
      .causalgrn_stop_cluster(cl, pids, hard = TRUE)
      stop("Interrupted.", call. = FALSE)
    }
  )
}
