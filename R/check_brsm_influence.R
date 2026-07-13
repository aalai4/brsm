#' Influence and Outlier Diagnostics for a brsm Fit
#'
#' Computes observation-level outlier and influence diagnostics for Bayesian
#' response surface models.
#'
#' Outlier diagnostics combine two signals:
#' \itemize{
#'   \item standardized residual magnitude (\code{|r_i| / sd(r)})
#'   \item posterior predictive interval coverage
#' }
#'
#' Influence diagnostics are based on analytical leverage values
#' \eqn{h_i = \text{diag}(X \mathbf{V}_n X^{\top})} computed from the
#' posterior covariance. Observations with high leverage are flagged
#' as influential.
#'
#' @param object A \code{brsm_fit} object from [fit_brsm()] or a
#'   \code{brsm_conjugate_fit} object.
#' @param ndraws Number of posterior predictive draws used to compute
#'   predictive means and intervals.
#' @param probs Length-2 numeric vector of predictive interval probabilities.
#' @param outlier_sd_threshold Threshold for standardized residual outliers.
#' @param leverage_threshold Leverage threshold above which observations
#'   are flagged as influential. If \code{NULL}, defaults to \code{3 * p / n}.
#' @param seed Optional random seed.
#' @param include_plot Logical; if \code{TRUE}, includes a ggplot residual
#'   diagnostic plot in the output.
#' @param verbose Logical; if \code{TRUE}, prints a concise summary.
#' @param ... Additional arguments.
#'
#' @return A list with components:
#'   \code{overview} (one-row data frame), \code{observations}
#'   (observation-level diagnostics), \code{passed} (logical scalar), and
#'   optionally \code{plot}.
#'
#' @examples
#' \dontrun{
#' fit <- fit_brsm(dat, response = "y", factor_names = c("x1", "x2"))
#' diag <- check_brsm_influence(fit)
#' diag$overview
#' head(diag$observations)
#' }
#' @export
check_brsm_influence <- function(object,
                                 ndraws = 400,
                                 probs = c(0.025, 0.975),
                                 outlier_sd_threshold = 3,
                                 leverage_threshold = NULL,
                                 seed = NULL,
                                 include_plot = FALSE,
                                 verbose = TRUE,
                                 ...) {
  if (!is.numeric(ndraws) || length(ndraws) != 1L || !is.finite(ndraws) ||
      ndraws < 1) {
    stop("ndraws must be a positive finite integer.")
  }
  ndraws <- as.integer(ndraws)

  probs <- .brsm_validate_probs(probs, require_length = 2)

  if (!is.numeric(outlier_sd_threshold) || length(outlier_sd_threshold) != 1L ||
      !is.finite(outlier_sd_threshold) || outlier_sd_threshold <= 0) {
    stop("outlier_sd_threshold must be a positive finite numeric scalar.")
  }

  if (!is.null(leverage_threshold) &&
      (!is.numeric(leverage_threshold) || length(leverage_threshold) != 1L ||
       !is.finite(leverage_threshold) || leverage_threshold <= 0)) {
    stop("leverage_threshold must be a positive finite numeric scalar.")
  }

  fit <- object
  if (inherits(object, "brsm_fit")) {
    fit <- object$fit
  }

  if (!inherits(fit, "brsm_conjugate_fit")) {
    stop("object must be a brsm_fit or brsm_conjugate_fit object in check_brsm_influence().")
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  y <- fit$y
  yrep <- posterior_predict_brsm(fit, max_draws = ndraws)

  if (!is.matrix(yrep) && !is.data.frame(yrep)) {
    stop("posterior_predict_brsm did not return a matrix-like object.")
  }
  yrep <- as.matrix(yrep)

  if (ncol(yrep) != length(y)) {
    stop("posterior_predict_brsm output columns do not match number of observations.")
  }

  yhat <- colMeans(yrep)
  lower <- apply(yrep, 2, stats::quantile, probs = probs[1], na.rm = TRUE)
  upper <- apply(yrep, 2, stats::quantile, probs = probs[2], na.rm = TRUE)

  residual <- y - yhat
  residual_scale <- .brsm_residual_scale(residual)
  std_residual <- residual / residual_scale

  outlier_std <- abs(std_residual) > outlier_sd_threshold
  outlier_interval <- y < lower | y > upper
  outlier_any <- outlier_std | outlier_interval

  X <- fit$X
  V_n <- fit$V_n
  leverage <- rowSums((X %*% V_n) * X)

  if (is.null(leverage_threshold)) {
    leverage_threshold <- 3 * ncol(X) / nrow(X)
  }
  influential <- leverage > leverage_threshold

  obs_diag <- data.frame(
    obs_id = seq_along(y),
    observed = y,
    predicted_mean = yhat,
    pred_lower = lower,
    pred_upper = upper,
    residual = residual,
    std_residual = std_residual,
    outlier_std = outlier_std,
    outlier_interval = outlier_interval,
    outlier = outlier_any,
    leverage = leverage,
    influential = influential,
    stringsAsFactors = FALSE
  )

  n_influential <- sum(obs_diag$influential, na.rm = TRUE)
  n_outliers <- sum(obs_diag$outlier, na.rm = TRUE)

  passed_outlier <- n_outliers == 0
  passed_influence <- n_influential == 0
  passed <- passed_outlier && passed_influence

  overview <- data.frame(
    n_obs = length(y),
    ndraws = ndraws,
    outlier_sd_threshold = outlier_sd_threshold,
    leverage_threshold = leverage_threshold,
    interval_lower_prob = probs[1],
    interval_upper_prob = probs[2],
    residual_sd = residual_scale,
    max_abs_std_residual = max(abs(std_residual), na.rm = TRUE),
    n_outliers = n_outliers,
    n_outlier_std = sum(obs_diag$outlier_std, na.rm = TRUE),
    n_outlier_interval = sum(obs_diag$outlier_interval, na.rm = TRUE),
    max_leverage = max(obs_diag$leverage, na.rm = TRUE),
    n_influential = n_influential,
    passed = passed,
    stringsAsFactors = FALSE
  )

  if (isTRUE(verbose)) {
    message(
      "check_brsm_influence: passed=", passed,
      ", n_outliers=", n_outliers,
      ", max_abs_std_resid=", format(overview$max_abs_std_residual, digits = 4),
      ", n_influential=", n_influential,
      ", max_leverage=", format(overview$max_leverage, digits = 4)
    )
  }

  out <- list(
    overview = overview,
    observations = obs_diag,
    passed = passed
  )

  if (isTRUE(include_plot)) {
    out$plot <- .brsm_plot_influence_diagnostics(
      obs_diag = obs_diag,
      outlier_sd_threshold = outlier_sd_threshold,
      leverage_threshold = leverage_threshold
    )
  }

  out
}


.brsm_residual_scale <- function(residuals) {
  if (!is.numeric(residuals) || length(residuals) == 0L) {
    return(1)
  }

  s <- stats::sd(residuals, na.rm = TRUE)
  if (is.finite(s) && s > 0) {
    return(s)
  }

  mad_scale <- stats::mad(residuals, center = stats::median(residuals, na.rm = TRUE),
    constant = 1.4826, na.rm = TRUE)
  if (is.finite(mad_scale) && mad_scale > 0) {
    return(mad_scale)
  }

  1
}


.brsm_plot_influence_diagnostics <- function(obs_diag,
                                             outlier_sd_threshold,
                                             leverage_threshold) {
  category <- ifelse(
    obs_diag$influential & obs_diag$outlier,
    "outlier+influential",
    ifelse(
      obs_diag$influential,
      "influential",
      ifelse(obs_diag$outlier, "outlier", "ok")
    )
  )

  plot_df <- obs_diag
  plot_df$category <- factor(
    category,
    levels = c("ok", "outlier", "influential", "outlier+influential")
  )

  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = predicted_mean, y = std_residual, color = category)
  ) +
    ggplot2::geom_point(alpha = 0.8, size = 2) +
    ggplot2::geom_hline(
      yintercept = c(-outlier_sd_threshold, outlier_sd_threshold),
      linetype = "dashed"
    ) +
    ggplot2::labs(
      title = "Influence and Outlier Diagnostics",
      subtitle = paste0("Leverage threshold = ", signif(leverage_threshold, 3)),
      x = "Posterior predictive mean",
      y = "Standardized residual",
      color = "Diagnostic"
    ) +
    ggplot2::theme_minimal()
}