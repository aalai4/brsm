# Centralized stationary-point solver status definitions.
.brsm_stationary_status_levels <- function() {
  c(
    "ok",
    "lapack_fail",
    "invalid_lu_diag",
    "kappa_exceeded"
  )
}

.brsm_stationary_status_labels <- function(status_code) {
  levels <- .brsm_stationary_status_levels()
  idx <- as.integer(status_code)
  labels <- rep(NA_character_, length(idx))
  valid <- !is.na(idx) & idx >= 0L & idx < length(levels)
  labels[valid] <- levels[idx[valid] + 1L]
  labels
}

.brsm_stationary_status_counts <- function(status_code) {
  levels <- .brsm_stationary_status_levels()
  idx <- as.integer(status_code)
  valid <- !is.na(idx) & idx >= 0L & idx < length(levels)
  stats::setNames(
    as.integer(tabulate(idx[valid] + 1L, nbins = length(levels))),
    levels
  )
}

#' Decode Stationary Solver Status Codes
#'
#' Converts per-draw stationary-point solver status codes into an ordered factor
#' with descriptive labels.
#'
#' Status meanings are:
#' \itemize{
#'   \item \code{0} = \code{"ok"}
#'   \item \code{1} = \code{"lapack_fail"}
#'   \item \code{2} = \code{"invalid_lu_diag"}
#'   \item \code{3} = \code{"kappa_exceeded"}
#' }
#'
#' @param status_code Integer-like vector of stationary solver status codes.
#'
#' @return An ordered factor with levels \code{"ok"}, \code{"lapack_fail"},
#'   \code{"invalid_lu_diag"}, and \code{"kappa_exceeded"}. Invalid or
#'   missing codes are returned as \code{NA}.
#'
#' @export
decode_stationary_status <- function(status_code) {
  labels <- .brsm_stationary_status_labels(status_code)
  factor(
    labels,
    levels = .brsm_stationary_status_levels(),
    ordered = TRUE
  )
}