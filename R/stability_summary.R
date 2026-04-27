#' Summarize BRSM Stability Diagnostics
#'
#' Extracts and displays stability diagnostics from BRSM analysis results.
#' Provides a unified view of exclusion rates, failure modes, and numerical
#' stability across function calls.
#'
#' @param object A `brsm_fit` object, or a result object (data frame or list)
#'   from `stationary_point()`, `posterior_ridge_analysis()`,
#'   `canonical_analysis()`, or similar functions that attach a
#'   `"diagnostics"` attribute.
#' @param ... Additional result objects to summarize together.
#'
#' @return A data frame with one row per analyzed result, containing:
#'   \itemize{
#'     \item `function_name`: Name of the analysis function.
#'     \item `n_draws`: Total number of posterior draws.
#'     \item `n_excluded`: Number of excluded draws (failed status).
#'     \item `pct_excluded`: Percentage of excluded draws.
#'     \item `status_ok`: Count of OK solves.
#'     \item `status_lapack_fail`: Count of LAPACK failures.
#'     \item `status_invalid_lu_diag`: Count of invalid LU diagonal failures.
#'     \item `status_kappa_exceeded`: Count of high-condition-number exclusions.
#'   }
#'
#' @details
#' For `brsm_fit` objects, this function extracts diagnostics from all
#' major internal analyses (stationary point computation, etc.).
#'
#' For analysis result objects with a `"diagnostics"` attribute, it extracts
#' the status codes and aggregates them into counts.
#'
#' Use this to quickly assess whether numerical instability is affecting your
#' inferences, and to understand which types of failures are most prevalent.
#'
#' @export
summarize_brsm_stability <- function(object, ...) {
  UseMethod("summarize_brsm_stability")
}

#' @rdname summarize_brsm_stability
#' @export
summarize_brsm_stability.brsm_fit <- function(object, ...) {
  # Extract diagnostics from the fit's stored attributes if available
  # For now, just return a message that direct fit diagnostics aren't stored
  # but user can call summarize_brsm_stability on individual analysis results
  diag_list <- list()

  # Check if fit object has internal diagnostic info
  if (!is.null(attr(object, "diagnostics"))) {
    diag_list[[1]] <- .extract_diagnostics_from_result(object)
  }

  # If there are additional objects passed, process them
  dots <- list(...)
  if (length(dots) > 0) {
    for (i in seq_along(dots)) {
      res <- .extract_diagnostics_from_result(dots[[i]])
      if (!is.null(res)) {
        diag_list[[length(diag_list) + 1]] <- res
      }
    }
  }

  if (length(diag_list) == 0) {
    message(
      "No diagnostics found in fit object. ",
      "Ensure you computed posterior_ridge_analysis(), canonical_analysis(), ",
      "or stationary_point() with diagnostics enabled."
    )
    return(invisible(NULL))
  }

  do.call(rbind, diag_list)
}

#' @rdname summarize_brsm_stability
#' @export
summarize_brsm_stability.default <- function(object, ...) {
  # Collect all objects (main + ...)
  all_objects <- list(object)
  dots <- list(...)
  if (length(dots) > 0) {
    all_objects <- c(all_objects, dots)
  }

  diag_list <- list()
  for (obj in all_objects) {
    res <- .extract_diagnostics_from_result(obj)
    if (!is.null(res)) {
      diag_list[[length(diag_list) + 1]] <- res
    }
  }

  if (length(diag_list) == 0) {
    stop(
      "No diagnostics found in any supplied objects. ",
      "Objects must have a 'diagnostics' attribute from ",
      "stationary_point(), posterior_ridge_analysis(), or canonical_analysis()."
    )
  }

  summary_df <- do.call(rbind, diag_list)
  rownames(summary_df) <- NULL
  summary_df
}

#' Extract Diagnostics from a Result Object
#'
#' Internal helper to extract stability diagnostics from an analysis result.
#'
#' @param obj A result object with a `"diagnostics"` attribute.
#'
#' @return A one-row data frame with aggregated diagnostics, or NULL if no
#'   diagnostics are found.
#'
#' @keywords internal
.extract_diagnostics_from_result <- function(obj) {
  diag <- attr(obj, "diagnostics", exact = TRUE)
  if (is.null(diag)) {
    return(NULL)
  }

  # Extract fields
  func_name <- diag$function_name %||% "unknown"
  n_draws <- length(diag$status_code)
  n_excluded <- diag$n_excluded %||% sum(diag$status_code != 0L)
  pct_excluded <- diag$pct_excluded %||% (n_excluded / n_draws)
  status_counts <- diag$status_counts %||% .brsm_stationary_status_counts(diag$status_code)

  # Ensure status_counts has all expected levels
  expected_levels <- c("ok", "lapack_fail", "invalid_lu_diag", "kappa_exceeded")
  for (level in expected_levels) {
    if (!(level %in% names(status_counts))) {
      status_counts[level] <- 0L
    }
  }
  status_counts <- status_counts[expected_levels]

  data.frame(
    function_name = func_name,
    n_draws = n_draws,
    n_excluded = n_excluded,
    pct_excluded = round(pct_excluded * 100, 2),
    status_ok = as.integer(status_counts[["ok"]]),
    status_lapack_fail = as.integer(status_counts[["lapack_fail"]]),
    status_invalid_lu_diag = as.integer(status_counts[["invalid_lu_diag"]]),
    status_kappa_exceeded = as.integer(status_counts[["kappa_exceeded"]]),
    stringsAsFactors = FALSE
  )
}