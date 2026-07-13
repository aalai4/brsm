#' Print Method for brsm_fit Objects
#'
#' Display a concise summary of a fitted Bayesian response surface model.
#'
#' @param x An object of class \code{brsm_fit}.
#' @param ... Additional arguments (unused).
#'
#' @return Invisibly returns \code{x}.
#'
#' @keywords internal
#' @export
print.brsm_fit <- function(x, ...) {
  cat("\nBayesian Response Surface Model\n\n")

  # Model specification
  cat("Model Formula:\n")
  cat("  ", format(x$formula), "\n\n")

  # Response and factors
  cat("Response Variable: ", x$response, "\n")
  cat("Factor Variables:  ", paste(x$factor_names, collapse = ", "), "\n")

  if (!is.null(x$coding) && is.list(x$coding) && !is.null(x$coding$method)) {
    cat("Coding Method:     ", x$coding$method, "\n")
    cat("Coding Stored:     yes\n")
  }

  cat("\n")

  # Factor ranges
  if (!is.null(x$ranges)) {
    cat("Factor Ranges:\n")
    for (fname in x$factor_names) {
      if (fname %in% names(x$ranges)) {
        rng <- x$ranges[[fname]]
        cat("  ", fname, ": [", rng[1], ", ", rng[2], "]\n", sep = "")
      }
    }
    cat("\n")
  }

  # Sampling information
  if (!is.null(x$sampling)) {
    cat("Conjugate Draws:\n")
    cat("  Draws:                  ", x$sampling$draws, "\n")
    cat("\n")
  }

  # Brief fit diagnostics
  if (!is.null(x$fit)) {
    cat("Model Fit Summary:\n")
    cat("  Observations:           ", nrow(x$fit$data), "\n")
    n_coef <- ncol(x$fit$draws) - 1
    cat("  Parameters estimated:   ", n_coef, "\n")

    ratio_message <- .brsm_low_information_message(
      n_obs = nrow(x$fit$data),
      n_coef = n_coef
    )
    if (!is.null(ratio_message)) {
      cat("  ", ratio_message, "\n", sep = "")
    }
    cat("\n")
  }

  cat("Call:\n")
  cat("  ")
  print(x$call)
  cat("\n")

  invisible(x)
}


#' Summary Method for brsm_fit Objects
#'
#' Display a detailed summary of a fitted Bayesian response surface model,
#' including coefficient summaries from the exact conjugate fit.
#'
#' @param object An object of class \code{brsm_fit}.
#' @param ... Additional arguments.
#'
#' @return An object of class \code{summary.brsm_fit} containing:
#'   \code{brsm_fit_obj} (the original brsm_fit object) and
#'   \code{fit_summary} (summary of the underlying conjugate fit).
#'
#' @keywords internal
#' @export
summary.brsm_fit <- function(object, ...) {
  # Extract and compute summary of underlying conjugate fit
  fit_obj <- object$fit
  fit_summary <- summary(fit_obj, ...)
  coef_uncertainty <- .brsm_fixed_effect_uncertainty(fit_obj)

  # Build custom summary object
  result <- list(
    brsm_fit_obj = object,
    fit_summary = fit_summary,
    coef_uncertainty = coef_uncertainty,
    formula = object$formula,
    response = object$response,
    factor_names = object$factor_names,
    coding = object$coding,
    ranges = object$ranges,
    sampling = object$sampling
  )

  class(result) <- c("summary.brsm_fit", "list")
  result
}


#' Print Method for summary.brsm_fit Objects
#'
#' Display a formatted summary of a Bayesian response surface model.
#'
#' @param x An object of class \code{summary.brsm_fit}.
#' @param ... Additional arguments (unused).
#'
#' @return Invisibly returns \code{x}.
#'
#' @keywords internal
#' @export
print.summary.brsm_fit <- function(x, ...) {
  cat("\nBayesian Response Surface Model Summary\n\n")

  # Model specification
  cat("Model Formula:\n")
  cat("  ", format(x$formula), "\n\n")

  cat("Response:     ", x$response, "\n")
  cat("Factors:      ", paste(x$factor_names, collapse = ", "), "\n")

  if (!is.null(x$coding) && is.list(x$coding) && !is.null(x$coding$method)) {
    cat("Coding:       ", x$coding$method, " (stored)\n")
  }

  cat("\n")

  # Factor ranges
  if (!is.null(x$ranges)) {
    cat("Factor Ranges:\n")
    for (fname in x$factor_names) {
      if (fname %in% names(x$ranges)) {
        rng <- x$ranges[[fname]]
        cat("  ", fname, ": [", rng[1], ", ", rng[2], "]\n", sep = "")
      }
    }
    cat("\n")
  }

  # Draws information
  if (!is.null(x$sampling)) {
    cat("Conjugate Draws:\n")
    cat("  Draws:                  ", x$sampling$draws, "\n")
    cat("\n")
  }

  if (!is.null(x$coef_uncertainty) && nrow(x$coef_uncertainty) > 0L) {
    n_overlap <- sum(x$coef_uncertainty$overlap_zero, na.rm = TRUE)
    median_pd <- stats::median(x$coef_uncertainty$pd, na.rm = TRUE)

    cat("Coefficient Uncertainty Check:\n")
    cat("  Median probability of direction: ",
      format(round(median_pd, 3), nsmall = 3), "\n",
      sep = ""
    )
    cat("  ", n_overlap, " of ", nrow(x$coef_uncertainty),
      " coefficients have 95% intervals overlapping 0.\n",
      sep = ""
    )
    if (n_overlap > 0L) {
      cat("  Caution: rank-ordering small coefficients is unstable when intervals overlap 0.\n")
    }
    cat("\n")
  }

  # Delegate to summary printing
  cat("Coefficient Summary:\n\n")
  print(x$fit_summary)

  invisible(x)
}