#' Start a cluster whose workers run the same SCMs code as this session
#'
#' Workers started as fresh R processes (PSOCK) load the *installed* SCMs, which differs from
#' the session's code whenever the package was loaded with `devtools::load_all()` or the
#' installed copy is stale. On Unix a fork cluster inherits the session's loaded namespace, so
#' the code always matches. Elsewhere, fail hard if the workers would run different code.
#'
#' @param n_cores Integer. Number of workers.
#' @param type Character. `"FORK"` (default on Unix) or `"PSOCK"`.
#' @return A cluster object; the caller stops it.
#' @noRd
make_scms_cluster <- function(n_cores, type = if (.Platform$OS.type == "unix") "FORK" else "PSOCK") {
  if (!type %in% c("FORK", "PSOCK")) {
    stop("type must be 'FORK' or 'PSOCK'.")
  }
  if (type == "FORK") {
    return(parallel::makeForkCluster(n_cores))
  }

  cl <- parallel::makeCluster(n_cores)
  session_path <- normalizePath(getNamespaceInfo(asNamespace("SCMs"), "path"))
  worker_paths <- tryCatch(
    unique(unlist(parallel::clusterEvalQ(
      cl, normalizePath(getNamespaceInfo(asNamespace("SCMs"), "path"))
    ))),
    error = function(e) {
      parallel::stopCluster(cl)
      stop("Parallel workers could not load SCMs: ", conditionMessage(e),
           ". Install the package or use cores = 1.", call. = FALSE)
    }
  )
  if (!identical(worker_paths, session_path)) {
    parallel::stopCluster(cl)
    stop("Parallel workers would run SCMs from ", paste(worker_paths, collapse = ", "),
         " but this session runs SCMs from ", session_path,
         ". Install this version of the package or use cores = 1.")
  }
  cl
}
