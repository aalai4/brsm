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

  # MCMC and sampling information
  if (!is.null(x$sampling)) {
    cat("MCMC Sampling:\n")
    cat("  Chains:                 ", x$sampling$chains, "\n")
    cat("  Iterations per chain:   ", x$sampling$iter, "\n")
    cat("  Warmup iterations:      ", x$sampling$warmup, "\n")
    cat("  Sampling preset:        ", x$sampling$sampling_preset, "\n")
    cat("\n")
  }

  # Brief fit diagnostics
  if (!is.null(x$fit)) {
    cat("Model Fit Summary:\n")
    cat("  Observations:           ", nrow(x$fit$data), "\n")
    cat("  Parameters estimated:   ", length(brms::fixef(x$fit)[, 1]), "\n")

    ratio_message <- .brsm_low_information_message(
      n_obs = nrow(x$fit$data),
      n_coef = length(brms::fixef(x$fit)[, 1])
    )
    if (!is.null(ratio_message)) {
      cat("  ", ratio_message, "\n", sep = "")
    }

    # Try to extract Rhat summary
    tryCatch(
      {
        rhat_vals <- brms::rhat(x$fit)
        if (!is.null(rhat_vals) && length(rhat_vals) > 0) {
          rhat_max <- max(rhat_vals, na.rm = TRUE)
          cat("  Max Rhat (convergence):  ", round(rhat_max, 4), "\n")
        }
      },
      error = function(e) NULL
    )
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
#' including MCMC diagnostics and coefficient summaries from the underlying
#' \code{brmsfit} object.
#'
#' @param object An object of class \code{brsm_fit}.
#' @param ... Additional arguments passed to \code{brms::summary.brmsfit()}.
#'
#' @return An object of class \code{summary.brsm_fit} containing:
#'   \code{brsm_fit_obj} (the original brsm_fit object) and
#'   \code{brmsfit_summary} (summary of the underlying brmsfit).
#'
#' @keywords internal
#' @export
summary.brsm_fit <- function(object, ...) {
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("package 'brms' is required for summary.brsm_fit().")
  }

  # Extract and compute summary of underlying brmsfit
  brmsfit_obj <- object$fit
  brmsfit_summary <- summary(brmsfit_obj, ...)
  coef_uncertainty <- .brsm_fixed_effect_uncertainty(brmsfit_obj)

  # Build custom summary object
  result <- list(
    brsm_fit_obj = object,
    brmsfit_summary = brmsfit_summary,
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

  # MCMC information
  if (!is.null(x$sampling)) {
    cat("MCMC Configuration:\n")
    cat("  ", x$sampling$chains, " chains x ", x$sampling$iter,
      " iterations (", x$sampling$warmup, " warmup)\n",
      sep = ""
    )
    cat("  Sampling profile: ", x$sampling$sampling_preset, "\n")
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

  # Delegate to brms summary printing
  cat("Coefficient Summary (via brms):\n\n")
  print(x$brmsfit_summary)

  invisible(x)
}