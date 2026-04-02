#' Posterior Predictive Check Summary for a brsm Fit
#'
#' Computes a compact posterior predictive summary using simulated responses
#' from the fitted model.
#'
#' @param object A \code{brsm_fit} object from [fit_brsm()] or a
#'   \code{brmsfit} object.
#' @param ndraws Number of posterior predictive draws.
#' @param probs Length-2 numeric vector of predictive interval probabilities.
#' @param seed Optional random seed.
#' @param include_plot Logical; if \code{TRUE}, includes a \code{pp_check}
#'   histogram overlay plot in the output.
#' @param ... Additional arguments passed to \code{brms::posterior_predict()}.
#'
#' @return A list with components:
#'   \code{summary} (one-row data frame), \code{observed} (numeric vector),
#'   \code{predicted_mean} (posterior mean prediction per observation), and
#'   optionally \code{plot}.
#'
#' @examples
#' \dontrun{
#' fit <- fit_brsm(dat, response = "y", factor_names = c("x1", "x2"))
#' ppc <- check_brsm_ppc(fit, ndraws = 200)
#' ppc$summary
#' }
#' @export
check_brsm_ppc <- function(object,
                           ndraws = 200,
                           probs = c(0.025, 0.975),
                           seed = NULL,
                           include_plot = FALSE,
                           ...) {
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("package 'brms' is required for check_brsm_ppc().")
  }

  fit <- .brsm_extract_fit(object, caller = "check_brsm_ppc")

  if (!is.numeric(ndraws) || length(ndraws) != 1 ||
        ndraws < 1 || !is.finite(ndraws)) {
    stop("ndraws must be a positive finite integer.")
  }
  ndraws <- as.integer(ndraws)

  probs <- .brsm_validate_probs(probs, require_length = 2)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  y <- as.numeric(brms::get_y(fit))
  yrep <- brms::posterior_predict(fit, ndraws = ndraws, ...)

  if (!is.matrix(yrep) && !is.data.frame(yrep)) {
    stop("posterior_predict did not return a matrix-like object.")
  }
  yrep <- as.matrix(yrep)

  yhat <- colMeans(yrep)
  lower <- apply(yrep, 2, stats::quantile, probs = probs[1], na.rm = TRUE)
  upper <- apply(yrep, 2, stats::quantile, probs = probs[2], na.rm = TRUE)
  coverage <- mean(y >= lower & y <= upper, na.rm = TRUE)

  rmse <- sqrt(mean((y - yhat)^2, na.rm = TRUE))

  summary_df <- data.frame(
    n_obs = length(y),
    ndraws = ndraws,
    observed_mean = mean(y, na.rm = TRUE),
    predicted_mean = mean(yhat, na.rm = TRUE),
    observed_sd = stats::sd(y, na.rm = TRUE),
    predicted_sd = stats::sd(yhat, na.rm = TRUE),
    interval_lower_prob = probs[1],
    interval_upper_prob = probs[2],
    interval_coverage = coverage,
    rmse = rmse,
    stringsAsFactors = FALSE
  )

  out <- list(
    summary = summary_df,
    observed = y,
    predicted_mean = yhat
  )

  if (isTRUE(include_plot)) {
    out$plot <- brms::pp_check(fit, ndraws = min(ndraws, 100), type = "hist")
  }

  out
}