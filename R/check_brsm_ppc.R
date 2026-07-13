#' Posterior Predictive Check Summary for a brsm Fit
#'
#' Computes a compact posterior predictive summary using simulated responses
#' from the fitted model.
#'
#' @param object A \code{brsm_fit} object from [fit_brsm()] or a
#'   \code{brsm_conjugate_fit} object.
#' @param ndraws Number of posterior predictive draws.
#' @param probs Length-2 numeric vector of predictive interval probabilities.
#' @param seed Optional random seed.
#' @param include_plot Logical; if \code{TRUE}, includes a ggplot2 density overlay.
#' @param ... Additional arguments.
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
  fit <- object
  if (inherits(object, "brsm_fit")) {
    fit <- object$fit
  }

  if (!inherits(fit, "brsm_conjugate_fit")) {
    stop("object must be a brsm_fit or brsm_conjugate_fit object.")
  }

  if (!is.numeric(ndraws) || length(ndraws) != 1 ||
        ndraws < 1 || !is.finite(ndraws)) {
    stop("ndraws must be a positive finite integer.")
  }
  ndraws <- as.integer(ndraws)

  probs <- .brsm_validate_probs(probs, require_length = 2)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  y <- fit$y
  yrep <- posterior_predict_brsm(fit, max_draws = ndraws)

  if (!is.matrix(yrep) && !is.data.frame(yrep)) {
    stop("posterior_predict_brsm did not return a matrix-like object.")
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
    n_plots <- min(ndraws, 50)
    plot_draws <- yrep[seq_len(n_plots), , drop = FALSE]
    
    # Reshape for ggplot
    df_y <- data.frame(value = y, group = "Observed", draw = 0L)
    df_yrep_list <- lapply(seq_len(n_plots), function(i) {
      data.frame(value = plot_draws[i, ], group = "Simulated", draw = i)
    })
    df_plot <- do.call(rbind, c(list(df_y), df_yrep_list))
    df_plot$group <- factor(df_plot$group, levels = c("Simulated", "Observed"))
    
    out$plot <- ggplot2::ggplot(df_plot, ggplot2::aes(x = value)) +
      ggplot2::geom_density(
        ggplot2::aes(group = interaction(group, draw), color = group, size = group, alpha = group)
      ) +
      ggplot2::scale_color_manual(values = c("Observed" = "black", "Simulated" = "lightblue")) +
      ggplot2::scale_size_manual(values = c("Observed" = 1.2, "Simulated" = 0.5)) +
      ggplot2::scale_alpha_manual(values = c("Observed" = 1.0, "Simulated" = 0.4)) +
      ggplot2::labs(
        title = "Posterior Predictive Check",
        subtitle = paste0("Observed vs. ", n_plots, " simulated datasets"),
        x = "y",
        y = "Density",
        color = "", size = "", alpha = ""
      ) +
      ggplot2::theme_minimal()
  }

  out
}