#' Build Prior Specifications for BRSM Models
#'
#' Builds \code{brms} priors for quadratic response-surface models used by
#' [fit_brsm()]. Priors are generated for intercept, linear, interaction,
#' quadratic, and residual sigma terms according to
#' \code{model_terms}.
#'
#' @param factor_names Character vector of factor names.
#' @param model_terms Polynomial term specification. One of
#'   \code{"second_order"}, \code{"first_order"}, \code{"first_order_twi"},
#'   or \code{"pure_quadratic"}.
#' @param prior_profile Prior profile controlling default scale parameters.
#'   One of \code{"legacy"} (wide flat priors), \code{"regularized"}
#'   (moderately shrinking), or \code{"adaptive"} (data-scaled). User-supplied
#'   scale arguments always override profile defaults.
#' @param coefficient_family Prior family for intercept and slope terms.
#'   One of \code{"normal"} or \code{"student_t"}.
#' @param intercept_sd Base scale for intercept prior.
#' @param linear_sd Base scale for linear term priors.
#' @param interaction_sd Base scale for interaction term priors.
#' @param quadratic_sd Base scale for quadratic term priors.
#' @param sigma_scale Base scale for residual \code{sigma} prior.
#' @param student_df Degrees of freedom when
#'   \code{coefficient_family = "student_t"}.
#' @param include_intercept Logical; include intercept prior.
#' @param include_sigma Logical; include residual sigma prior.
#' @param autoscale Logical; if \code{TRUE} and \code{data}/\code{response}
#'   are provided, scales all prior standard deviations by
#'   \code{sd(data[[response]])}.
#' @param data Optional data frame used for autoscaling.
#' @param response Optional response column name used for autoscaling.
#'
#' @return A \code{brmsprior} object that can be passed directly to
#'   [fit_brsm()] or \code{brms::brm()} as the \code{prior} argument.
#' @export
specify_brsm_priors <- function(
    factor_names,
    model_terms = c(
      "second_order", "first_order", "first_order_twi", "pure_quadratic"
    ),
    prior_profile = c("legacy", "regularized", "adaptive"),
    coefficient_family = c("normal", "student_t"),
    intercept_sd = 5,
    linear_sd = 2,
    interaction_sd = 1,
    quadratic_sd = 1,
    sigma_scale = 2.5,
    student_df = 3,
    include_intercept = TRUE,
    include_sigma = TRUE,
    autoscale = FALSE,
    data = NULL,
    response = NULL) {
  factor_names <- .brsm_validate_factor_names(factor_names)
  model_terms <- match.arg(model_terms)
  prior_profile <- match.arg(prior_profile)
  coefficient_family <- match.arg(coefficient_family)

  intercept_sd_user <- !missing(intercept_sd)
  linear_sd_user <- !missing(linear_sd)
  interaction_sd_user <- !missing(interaction_sd)
  quadratic_sd_user <- !missing(quadratic_sd)
  sigma_scale_user <- !missing(sigma_scale)
  coefficient_family_user <- !missing(coefficient_family)
  student_df_user <- !missing(student_df)

  if (prior_profile != "legacy") {
    p <- length(factor_names)
    p_twi <- if (p > 1L) choose(p, 2L) else 0L

    # Profile defaults by polynomial order on coded predictors.
    prof <- switch(
      model_terms,
      first_order = list(
        intercept_sd = 3,
        linear_sd = 1.5,
        interaction_sd = 1,
        quadratic_sd = 1,
        sigma_scale = 2,
        coefficient_family = "student_t",
        student_df = 3,
        n_coef = p
      ),
      first_order_twi = list(
        intercept_sd = 3,
        linear_sd = 1.25,
        interaction_sd = 0.6,
        quadratic_sd = 1,
        sigma_scale = 2,
        coefficient_family = "student_t",
        student_df = 3,
        n_coef = p + p_twi
      ),
      pure_quadratic = list(
        intercept_sd = 3,
        linear_sd = 1.25,
        interaction_sd = 1,
        quadratic_sd = 0.6,
        sigma_scale = 2,
        coefficient_family = "student_t",
        student_df = 3,
        n_coef = 2L * p
      ),
      second_order = list(
        intercept_sd = 3,
        linear_sd = 1,
        interaction_sd = 0.5,
        quadratic_sd = 0.5,
        sigma_scale = 2,
        coefficient_family = "student_t",
        student_df = 3,
        n_coef = 2L * p + p_twi
      )
    )

    if (!intercept_sd_user) intercept_sd <- prof$intercept_sd
    if (!linear_sd_user) linear_sd <- prof$linear_sd
    if (!interaction_sd_user) interaction_sd <- prof$interaction_sd
    if (!quadratic_sd_user) quadratic_sd <- prof$quadratic_sd
    if (!sigma_scale_user) sigma_scale <- prof$sigma_scale
    if (!coefficient_family_user) coefficient_family <- prof$coefficient_family
    if (!student_df_user) student_df <- prof$student_df

    # Adaptive profile tightens slope priors when data are weak for the model.
    if (prior_profile == "adaptive" && !is.null(data) && is.data.frame(data)) {
      n <- nrow(data)
      info_ratio <- n / max(1L, prof$n_coef)

      if (is.finite(info_ratio) && info_ratio < 3) {
        linear_sd <- linear_sd * 0.8
        interaction_sd <- interaction_sd * 0.6
        quadratic_sd <- quadratic_sd * 0.6
      } else if (is.finite(info_ratio) && info_ratio < 5) {
        linear_sd <- linear_sd * 0.9
        interaction_sd <- interaction_sd * 0.75
        quadratic_sd <- quadratic_sd * 0.75
      }
    }
  }

  .brsm_check_positive_scalar(intercept_sd, "intercept_sd")
  .brsm_check_positive_scalar(linear_sd, "linear_sd")
  .brsm_check_positive_scalar(interaction_sd, "interaction_sd")
  .brsm_check_positive_scalar(quadratic_sd, "quadratic_sd")
  .brsm_check_positive_scalar(sigma_scale, "sigma_scale")
  .brsm_check_positive_scalar(student_df, "student_df")

  if (!is.logical(include_intercept) || length(include_intercept) != 1L ||
      is.na(include_intercept)) {
    stop("include_intercept must be TRUE or FALSE.")
  }
  if (!is.logical(include_sigma) || length(include_sigma) != 1L ||
      is.na(include_sigma)) {
    stop("include_sigma must be TRUE or FALSE.")
  }
  if (!is.logical(autoscale) || length(autoscale) != 1L || is.na(autoscale)) {
    stop("autoscale must be TRUE or FALSE.")
  }

  scale_multiplier <- 1
  if (isTRUE(autoscale)) {
    if (is.null(data) || !is.data.frame(data)) {
      stop("When autoscale=TRUE, data must be provided as a data.frame.")
    }
    if (is.null(response) || !is.character(response) || length(response) != 1L) {
      stop("When autoscale=TRUE, response must be a single character name.")
    }
    if (!response %in% names(data)) {
      stop("response column not found in data.")
    }
    if (!is.numeric(data[[response]])) {
      stop("response column must be numeric for autoscaling.")
    }

    y_sd <- stats::sd(data[[response]], na.rm = TRUE)
    if (!is.finite(y_sd) || y_sd <= 0) {
      stop("response standard deviation must be finite and > 0 for autoscaling.")
    }
    scale_multiplier <- y_sd
  }

  intercept_scale <- intercept_sd * scale_multiplier
  linear_scale <- linear_sd * scale_multiplier
  interaction_scale <- interaction_sd * scale_multiplier
  quadratic_scale <- quadratic_sd * scale_multiplier
  sigma_prior_scale <- sigma_scale * scale_multiplier

  linear_terms <- factor_names
  interaction_terms <- character(0)
  if (length(factor_names) > 1L) {
    pairs <- utils::combn(factor_names, 2, simplify = FALSE)
    interaction_terms <- vapply(pairs, function(pair) {
      paste0(pair[[1]], ":", pair[[2]])
    }, character(1))
  }
  quadratic_terms <- paste0("I(", factor_names, "^2)")

  include_interactions <- model_terms %in% c("second_order", "first_order_twi")
  include_quadratic <- model_terms %in% c("second_order", "pure_quadratic")

  available_b_coefs <- .brsm_default_prior_b_coefs(
    response = response,
    linear_terms = linear_terms,
    interaction_terms = interaction_terms,
    quadratic_terms = quadratic_terms,
    include_interactions = include_interactions,
    include_quadratic = include_quadratic,
    data = data
  )

  term_groups <- list(
    linear = linear_terms,
    interaction = if (include_interactions) interaction_terms else character(0),
    quadratic = if (include_quadratic) quadratic_terms else character(0)
  )
  resolved <- .brsm_resolve_b_prior_targets(term_groups, available_b_coefs)

  prior_list <- list()

  .brsm_set_prior <- function(prior_str, class_val, coef_val = "") {
    df <- data.frame(
      prior = prior_str,
      class = class_val,
      coef = coef_val,
      group = "",
      resp = "",
      dpar = "",
      nlpar = "",
      lb = NA_character_,
      ub = NA_character_,
      source = "user",
      stringsAsFactors = FALSE
    )
    class(df) <- c("brmsprior", "data.frame")
    df
  }

  if (isTRUE(include_intercept)) {
    prior_list[[length(prior_list) + 1L]] <- .brsm_set_prior(
      .brsm_prior_string(coefficient_family, 0, intercept_scale, student_df),
      class_val = "Intercept"
    )
  }

  need_global_b <- !is.null(available_b_coefs) && any(lengths(resolved$unmatched) > 0L)
  if (need_global_b) {
    prior_list[[length(prior_list) + 1L]] <- .brsm_set_prior(
      .brsm_prior_string(coefficient_family, 0, linear_scale, student_df),
      class_val = "b"
    )
  }

  for (f in resolved$matched$linear) {
    prior_list[[length(prior_list) + 1L]] <- .brsm_set_prior(
      .brsm_prior_string(coefficient_family, 0, linear_scale, student_df),
      class_val = "b",
      coef_val = f
    )
  }

  if (length(resolved$matched$interaction) > 0L) {
    for (coef_name in resolved$matched$interaction) {
      prior_list[[length(prior_list) + 1L]] <- .brsm_set_prior(
        .brsm_prior_string(coefficient_family, 0, interaction_scale, student_df),
        class_val = "b",
        coef_val = coef_name
      )
    }
  }

  if (length(resolved$matched$quadratic) > 0L) {
    for (coef_name in resolved$matched$quadratic) {
      prior_list[[length(prior_list) + 1L]] <- .brsm_set_prior(
        .brsm_prior_string(coefficient_family, 0, quadratic_scale, student_df),
        class_val = "b",
        coef_val = coef_name
      )
    }
  }

  if (isTRUE(include_sigma)) {
    prior_list[[length(prior_list) + 1L]] <- .brsm_set_prior(
      paste0("student_t(3, 0, ", signif(sigma_prior_scale, 6), ")"),
      class_val = "sigma"
    )
  }

  out <- do.call(rbind, prior_list)
  class(out) <- c("brmsprior", "data.frame")
  rownames(out) <- NULL
  out
}


.brsm_prior_string <- function(family, location, scale, student_df) {
  if (family == "normal") {
    return(paste0("normal(", location, ", ", signif(scale, 6), ")"))
  }

  paste0(
    "student_t(",
    signif(student_df, 6),
    ", ",
    location,
    ", ",
    signif(scale, 6),
    ")"
  )
}


.brsm_check_positive_scalar <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
    stop(name, " must be a finite numeric scalar > 0.")
  }
}


.brsm_default_prior_b_coefs <- function(response,
                                        linear_terms,
                                        interaction_terms,
                                        quadratic_terms,
                                        include_interactions,
                                        include_quadratic,
                                        data) {
  rhs_terms <- c(linear_terms)
  if (isTRUE(include_interactions)) {
    rhs_terms <- c(rhs_terms, interaction_terms)
  }
  if (isTRUE(include_quadratic)) {
    rhs_terms <- c(rhs_terms, quadratic_terms)
  }

  if (length(rhs_terms) == 0L) {
    return(NULL)
  }

  rhs_terms
}


.brsm_resolve_b_prior_targets <- function(term_groups, available_b_coefs) {
  if (is.null(available_b_coefs)) {
    return(list(
      matched = term_groups,
      unmatched = lapply(term_groups, function(x) character(0))
    ))
  }

  matched <- list(linear = character(0), interaction = character(0), quadratic = character(0))
  unmatched <- list(linear = character(0), interaction = character(0), quadratic = character(0))

  for (coef in term_groups$linear) {
    resolved <- .brsm_resolve_linear_coef(coef, available_b_coefs)
    if (is.na(resolved)) {
      unmatched$linear <- c(unmatched$linear, coef)
    } else {
      matched$linear <- c(matched$linear, resolved)
    }
  }

  for (coef in term_groups$interaction) {
    resolved <- .brsm_resolve_interaction_coef(coef, available_b_coefs)
    if (is.na(resolved)) {
      unmatched$interaction <- c(unmatched$interaction, coef)
    } else {
      matched$interaction <- c(matched$interaction, resolved)
    }
  }

  for (coef in term_groups$quadratic) {
    resolved <- .brsm_resolve_quadratic_coef(coef, available_b_coefs)
    if (is.na(resolved)) {
      unmatched$quadratic <- c(unmatched$quadratic, coef)
    } else {
      matched$quadratic <- c(matched$quadratic, resolved)
    }
  }

  matched <- lapply(matched, unique)
  list(matched = matched, unmatched = unmatched)
}


.brsm_resolve_linear_coef <- function(coef, available_b_coefs) {
  if (coef %in% available_b_coefs) {
    return(coef)
  }

  alt <- make.names(coef)
  if (alt %in% available_b_coefs) {
    return(alt)
  }

  NA_character_
}


.brsm_resolve_interaction_coef <- function(coef, available_b_coefs) {
  parts <- strsplit(coef, ":", fixed = TRUE)[[1]]
  if (length(parts) != 2L) {
    return(NA_character_)
  }

  c1 <- parts[1]
  c2 <- parts[2]
  candidates <- unique(c(
    paste0(c1, ":", c2),
    paste0(c2, ":", c1),
    paste0(c1, ".", c2),
    paste0(c2, ".", c1),
    make.names(paste0(c1, ":", c2)),
    make.names(paste0(c2, ":", c1))
  ))

  found <- intersect(candidates, available_b_coefs)
  if (length(found) == 0L) NA_character_ else found[[1L]]
}


.brsm_resolve_quadratic_coef <- function(coef, available_b_coefs) {
  f <- sub("^I\\((.+)\\^2\\)$", "\\1", coef)
  canonical <- paste0("I(", f, "^2)")
  candidates <- unique(c(
    canonical,
    make.names(canonical),
    paste0("I", f, "E2"),
    paste0("I", f, ".2"),
    paste0("I.", f, ".2")
  ))

  found <- intersect(candidates, available_b_coefs)
  if (length(found) > 0L) {
    return(found[[1L]])
  }

  idx <- which(
    grepl(f, available_b_coefs, fixed = TRUE) &
      (grepl("\\^2\\)", available_b_coefs) |
         grepl("E2", available_b_coefs, fixed = TRUE) |
         grepl("\\.2\\.?", available_b_coefs))
  )

  if (length(idx) == 0L) NA_character_ else available_b_coefs[idx[[1L]]]
}


#' Check Prior Specifications with Predictive Simulation
#'
#' Visualizes the prior predictive distribution implied by specified priors,
#' allowing users to assess whether priors place reasonable mass on plausible
#' response surfaces before committing to full Bayesian fitting.
#'
#' Samples coefficients and residual standard deviation from the prior
#' specification, then generates and visualizes predictions across factor
#' ranges. Warnings are raised if the prior predictive distribution seems
#' misaligned with observed response data.
#'
#' @param data Data frame containing response and factor columns.
#' @param response Name of the response variable (used for scale reference).
#' @param factor_names Character vector of factor names.
#' @param prior Optional prior specification. If \code{NULL}, uses
#'   \code{prior_profile} to generate priors.
#' @param prior_profile Prior profile to use when \code{prior = NULL}.
#'   One of \code{"legacy"}, \code{"regularized"}, or \code{"adaptive"}.
#'   Default is \code{"legacy"}.
#' @param model_terms Polynomial term specification (one of
#'   \code{"second_order"}, \code{"first_order"}, \code{"first_order_twi"},
#'   or \code{"pure_quadratic"}). Default is \code{"second_order"}.
#' @param n_prior_samples Number of prior samples to draw. Default is 100.
#' @param n_grid_per_factor Number of grid points per factor. Default is 10.
#' @param plot Logical; if \code{TRUE}, creates base-R plots of prior
#'   predictive distributions. Default is \code{TRUE}.
#' @param seed Optional random seed for reproducibility.
#'
#' @return Invisibly returns a list with elements:
#'   \code{prior_samples} (data frame of sampled coefficients and sigma),
#'   \code{predictions} (predictions across grid),
#'   \code{summary} (summary statistics including min/max/sd of predictions,
#'   and alignment flags).
#'   Also prints a summary table of predictions and any warnings.
#'
#' @examples
#' \dontrun{
#' dat <- data.frame(
#'   x1 = runif(50, -1, 1),
#'   x2 = runif(50, -1, 1),
#'   y = rnorm(50, mean = 5, sd = 1)
#' )
#' check_brsm_priors(
#'   data = dat,
#'   response = "y",
#'   factor_names = c("x1", "x2"),
#'   prior_profile = "regularized",
#'   n_prior_samples = 200,
#'   seed = 123
#' )
#' }
#'
#' @export
check_brsm_priors <- function(data,
                              response,
                              factor_names,
                              prior = NULL,
                              prior_profile = c("legacy", "regularized", "adaptive"),
                              model_terms = c("second_order", "first_order",
                                              "first_order_twi", "pure_quadratic"),
                              n_prior_samples = 100L,
                              n_grid_per_factor = 10L,
                              plot = TRUE,
                              seed = NULL) {
  if (!is.data.frame(data)) {
    stop("data must be a data.frame.")
  }
  if (!is.character(response) || length(response) != 1L) {
    stop("response must be a single character string.")
  }
  if (!response %in% names(data)) {
    stop("response '", response, "' not found in data.")
  }
  if (!is.numeric(data[[response]])) {
    stop("response column must be numeric.")
  }

  factor_names <- .brsm_validate_factor_names(factor_names)
  model_terms <- match.arg(model_terms)
  prior_profile <- match.arg(prior_profile)

  n_prior_samples <- as.integer(n_prior_samples)
  n_grid_per_factor <- as.integer(n_grid_per_factor)

  if (n_prior_samples < 1L) {
    stop("n_prior_samples must be >= 1.")
  }
  if (n_grid_per_factor < 2L) {
    stop("n_grid_per_factor must be >= 2.")
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Generate or validate prior
  if (is.null(prior)) {
    prior <- specify_brsm_priors(
      factor_names = factor_names,
      model_terms = model_terms,
      prior_profile = prior_profile,
      autoscale = TRUE,
      data = data,
      response = response
    )
  }

  # Sample from prior predictive
  prior_samples <- .brsm_sample_prior_predictive(
    prior = prior,
    factor_names = factor_names,
    model_terms = model_terms,
    n_samples = n_prior_samples
  )

  # Create prediction grid
  grid <- .brsm_create_prediction_grid(
    data = data,
    factor_names = factor_names,
    model_terms = model_terms,
    n_grid_per_factor = n_grid_per_factor
  )

  # Generate predictions
  predictions <- .brsm_predict_from_prior_samples(
    prior_samples = prior_samples,
    grid = grid,
    factor_names = factor_names,
    model_terms = model_terms
  )

  # Compute summary statistics
  y_obs <- data[[response]]
  y_obs_range <- range(y_obs, na.rm = TRUE)
  y_obs_mean <- mean(y_obs, na.rm = TRUE)
  y_obs_sd <- stats::sd(y_obs, na.rm = TRUE)

  pred_stats <- list(
    min = min(predictions$pred, na.rm = TRUE),
    max = max(predictions$pred, na.rm = TRUE),
    mean = mean(predictions$pred, na.rm = TRUE),
    sd = stats::sd(predictions$pred, na.rm = TRUE),
    q05 = stats::quantile(predictions$pred, 0.05, na.rm = TRUE),
    q95 = stats::quantile(predictions$pred, 0.95, na.rm = TRUE)
  )

  # Check for alignment issues
  warnings_list <- character(0)

  # Check 1: Does prior predictive range substantially exclude observed data?
  coverage_ratio <- abs(y_obs_range[2] - y_obs_range[1]) /
                     abs(pred_stats$q95 - pred_stats$q05)
  if (coverage_ratio > 2) {
    warnings_list <- c(warnings_list,
      "Prior predictive interval is much wider than observed range. Prior may be too diffuse.")
  }

  # Check 2: Does prior mean imply plausible baseline?
  if (abs(pred_stats$mean - y_obs_mean) > 2 * y_obs_sd) {
    warnings_list <- c(warnings_list,
      "Prior predictive mean differs substantially from observed mean. Consider adjusting intercept prior.")
  }

  # Create summary
  summary_table <- data.frame(
    Statistic = c("Observed Min", "Observed Max", "Observed Mean", "Observed SD",
                  "Prior Pred Min", "Prior Pred Max", "Prior Pred Mean", "Prior Pred SD",
                  "Prior Pred 5%", "Prior Pred 95%"),
    Value = c(
      y_obs_range[1], y_obs_range[2], y_obs_mean, y_obs_sd,
      pred_stats$min, pred_stats$max, pred_stats$mean, pred_stats$sd,
      pred_stats$q05, pred_stats$q95
    )
  )

  # Print summary
  cat("\nPrior Predictive Check:\n\n")
  cat("Model: ", model_terms, "\n")
  cat("Factors: ", paste(factor_names, collapse = ", "), "\n")
  cat("Prior samples drawn: ", n_prior_samples, "\n")
  cat("Prediction grid: ", paste(rep(n_grid_per_factor, length(factor_names)), collapse = " x "), "\n")
  cat("\n")
  print(summary_table, row.names = FALSE)
  cat("\n")

  if (length(warnings_list) > 0L) {
    cat("Warnings:\n")
    for (w in warnings_list) {
      cat("  - ", w, "\n", sep = "")
    }
    cat("\n")
  }

  # Visualization
  if (isTRUE(plot)) {
    .brsm_visualize_prior_predictive(
      predictions = predictions,
      factor_names = factor_names,
      y_obs = y_obs,
      model_terms = model_terms
    )
  }

  invisible(list(
    prior_samples = prior_samples,
    predictions = predictions,
    summary = list(
      stats = pred_stats,
      obs_range = y_obs_range,
      obs_mean = y_obs_mean,
      obs_sd = y_obs_sd,
      warnings = warnings_list
    )
  ))
}


.brsm_sample_prior_predictive <- function(prior, factor_names, model_terms, n_samples) {
  # Extract prior specifications
  prior_df <- as.data.frame(prior)

  # Determine coefficient names based on model_terms
  linear_coefs <- factor_names
  interaction_coefs <- if (length(factor_names) > 1) {
    pairs <- utils::combn(factor_names, 2, simplify = FALSE)
    vapply(pairs, function(p) paste0(p[1], ":", p[2]), character(1))
  } else character(0)
  quadratic_coefs <- if (model_terms %in% c("second_order", "pure_quadratic")) {
    paste0("I(", factor_names, "^2)")
  } else character(0)

  # Initialize sample matrix
  all_coefs <- c("Intercept", linear_coefs, interaction_coefs, quadratic_coefs)
  samples <- matrix(0, nrow = n_samples, ncol = length(all_coefs) + 1)
  colnames(samples) <- c(all_coefs, "sigma")

  # Sample from each prior distribution
  for (coef in all_coefs) {
    coef_rows <- which(prior_df$class == "b" & prior_df$coef == coef)
    if (length(coef_rows) == 0) {
      coef_rows <- which(prior_df$class == "Intercept" & coef == "Intercept")
    }

    if (length(coef_rows) == 0) {
      # Default: normal(0, 1) if no prior specified
      samples[, coef] <- stats::rnorm(n_samples, mean = 0, sd = 1)
    } else {
      prior_row <- prior_df[coef_rows[1], ]
      samples[, coef] <- .brsm_sample_from_prior_spec(prior_row, n_samples)
    }
  }

  # Sample sigma
  sigma_rows <- which(prior_df$class == "sigma")
  if (length(sigma_rows) > 0) {
    prior_row <- prior_df[sigma_rows[1], ]
    samples[, "sigma"] <- abs(.brsm_sample_from_prior_spec(prior_row, n_samples))
  } else {
    # Default: student_t(3, 0, 2.5)
    samples[, "sigma"] <- abs(stats::rt(n_samples, df = 3) * 2.5)
  }

  as.data.frame(samples)
}


# Sample from a single prior specification (internal helper)
.brsm_sample_from_prior_spec <- function(prior_row, n_samples) {
  prior_str <- prior_row$prior

  if (is.na(prior_str) || !nzchar(prior_str)) {
    return(stats::rnorm(n_samples, 0, 1))
  }

  # Parse common prior families
  # normal(mu, sigma)
  if (grepl("^normal\\(", prior_str)) {
    params <- .brsm_extract_prior_params(prior_str)
    mu <- if (length(params) > 0) params[1] else 0
    sigma <- if (length(params) > 1) params[2] else 1
    return(stats::rnorm(n_samples, mean = mu, sd = sigma))
  }

  # student_t(df, mu, sigma)
  if (grepl("^student_t\\(", prior_str)) {
    params <- .brsm_extract_prior_params(prior_str)
    df <- if (length(params) > 0) params[1] else 3
    mu <- if (length(params) > 1) params[2] else 0
    sigma <- if (length(params) > 2) params[3] else 1
    return(mu + sigma * stats::rt(n_samples, df = df))
  }

  # exponential(lambda)
  if (grepl("^exponential\\(", prior_str)) {
    params <- .brsm_extract_prior_params(prior_str)
    lambda <- if (length(params) > 0) params[1] else 1
    return(stats::rexp(n_samples, rate = lambda))
  }

  # Default fallback
  stats::rnorm(n_samples, 0, 1)
}


# Extract numeric parameters from a prior string (internal helper)
.brsm_extract_prior_params <- function(prior_str) {
  # Remove family name and parentheses
  inner <- sub("^[a-z_]+\\((.*)\\)$", "\\1", prior_str)
  if (inner == prior_str) {
    return(numeric(0))
  }

  # Split by comma and convert to numeric
  parts <- strsplit(inner, ",", fixed = TRUE)[[1]]
  as.numeric(trimws(parts))
}


# Create a prediction grid for factors
.brsm_create_prediction_grid <- function(data, factor_names, model_terms, n_grid_per_factor) {
  grids <- lapply(factor_names, function(f) {
    vals <- data[[f]]
    vals <- vals[!is.na(vals)]
    if (length(vals) == 0) {
      seq(-1, 1, length.out = n_grid_per_factor)
    } else {
      seq(min(vals), max(vals), length.out = n_grid_per_factor)
    }
  })
  names(grids) <- factor_names

  # Create full grid
  grid_df <- do.call(expand.grid, grids)

  # Compute derived terms
  linear_coefs <- factor_names

  interaction_terms <- if (model_terms %in% c("second_order", "first_order_twi") &&
                            length(factor_names) > 1) {
    pairs <- utils::combn(factor_names, 2, simplify = FALSE)
    for (pair in pairs) {
      col_name <- paste0(pair[1], ":", pair[2])
      grid_df[[col_name]] <- grid_df[[pair[1]]] * grid_df[[pair[2]]]
    }
    vapply(pairs, function(p) paste0(p[1], ":", p[2]), character(1))
  } else character(0)

  quadratic_terms <- if (model_terms %in% c("second_order", "pure_quadratic")) {
    for (f in factor_names) {
      col_name <- paste0("I(", f, "^2)")
      grid_df[[col_name]] <- grid_df[[f]]^2
    }
    paste0("I(", factor_names, "^2)")
  } else character(0)

  grid_df$grid_id <- seq_len(nrow(grid_df))
  grid_df
}


# Generate predictions from prior samples and grid
.brsm_predict_from_prior_samples <- function(prior_samples, grid, factor_names, model_terms) {
  linear_coefs <- factor_names
  interaction_coefs <- if (length(factor_names) > 1) {
    pairs <- utils::combn(factor_names, 2, simplify = FALSE)
    vapply(pairs, function(p) paste0(p[1], ":", p[2]), character(1))
  } else character(0)
  quadratic_coefs <- if (model_terms %in% c("second_order", "pure_quadratic")) {
    paste0("I(", factor_names, "^2)")
  } else character(0)

  all_coefs <- c(linear_coefs, interaction_coefs, quadratic_coefs)

  # Predictions: Intercept + sum of (coef * value) + normal error
  n_samples <- nrow(prior_samples)
  n_grid <- nrow(grid)

  pred_matrix <- matrix(0, nrow = n_grid, ncol = n_samples)

  for (i in seq_len(n_grid)) {
    # Intercept
    pred_matrix[i, ] <- prior_samples$Intercept

    # Linear terms
    for (coef in linear_coefs) {
      if (coef %in% names(grid)) {
        pred_matrix[i, ] <- pred_matrix[i, ] + prior_samples[[coef]] * grid[i, coef]
      }
    }

    # Interaction terms
    for (coef in interaction_coefs) {
      if (coef %in% names(grid)) {
        pred_matrix[i, ] <- pred_matrix[i, ] + prior_samples[[coef]] * grid[i, coef]
      }
    }

    # Quadratic terms
    for (coef in quadratic_coefs) {
      if (coef %in% names(grid)) {
        pred_matrix[i, ] <- pred_matrix[i, ] + prior_samples[[coef]] * grid[i, coef]
      }
    }

    # Add residual variability
    pred_matrix[i, ] <- pred_matrix[i, ] + stats::rnorm(n_samples, 0, prior_samples$sigma)
  }

  # Flatten predictions
  predictions_vec <- as.numeric(pred_matrix)

  # Replicate grid for each sample
  grid_rep <- grid[rep(seq_len(nrow(grid)), each = n_samples), ]
  sample_id <- rep(seq_len(n_samples), n_grid)

  data.frame(
    sample_id = sample_id,
    grid_rep,
    pred = predictions_vec
  )
}


# Visualize prior predictive distribution
.brsm_visualize_prior_predictive <- function(predictions, factor_names, y_obs, model_terms) {
  if (length(factor_names) == 1L) {
    # 1D case: plot posterior predictive density + observed histogram
    oldpar <- graphics::par(mfrow = c(1, 1))
    on.exit(graphics::par(oldpar))

    graphics::hist(y_obs, breaks = "Sturges", main = "Prior Predictive Check (1D)",
                   xlab = "Response", freq = FALSE, col = "lightgray", alpha = 0.7)
    graphics::lines(stats::density(predictions$pred, na.rm = TRUE), col = "blue", lwd = 2)
    graphics::legend("topright",
                     legend = c("Observed data", "Prior predictive"),
                     col = c("gray", "blue"),
                     lty = c(0, 1), pch = c(15, NA_integer_))

  } else if (length(factor_names) == 2L) {
    # 2D case: plot prior predictive surface for a few samples
    oldpar <- graphics::par(mfrow = c(2, 2))
    on.exit(graphics::par(oldpar))

    unique_samples <- unique(predictions$sample_id)[1:4]

    for (sample_idx in unique_samples) {
      pred_subset <- predictions[predictions$sample_id == sample_idx, ]
      x1_name <- factor_names[1]
      x2_name <- factor_names[2]

      x1_vals <- sort(unique(pred_subset[[x1_name]]))
      x2_vals <- sort(unique(pred_subset[[x2_name]]))

      z_mat <- matrix(pred_subset$pred, nrow = length(x1_vals), ncol = length(x2_vals))

      graphics::contour(x1_vals, x2_vals, z_mat,
                        main = paste("Prior Sample", sample_idx),
                        xlab = x1_name, ylab = x2_name)
    }

  } else {
    # High-dimensional case: just plot prior predictive density
    oldpar <- graphics::par(mfrow = c(1, 1))
    on.exit(graphics::par(oldpar))

    graphics::hist(y_obs, breaks = "Sturges",
                   main = paste("Prior Predictive Check (", length(factor_names), " factors)"),
                   xlab = "Response", freq = FALSE, col = "lightgray", alpha = 0.7)
    graphics::lines(stats::density(predictions$pred, na.rm = TRUE), col = "blue", lwd = 2)
    graphics::legend("topright",
                     legend = c("Observed data", "Prior predictive"),
                     col = c("gray", "blue"),
                     lty = c(0, 1), pch = c(15, NA_integer_))
  }

  invisible(NULL)
}