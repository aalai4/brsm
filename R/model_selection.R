#' Recommend Model Terms for Response Surface Models
#'
#' Analyzes data patterns and recommends appropriate polynomial terms for
#' BRSM model fitting. Fits exploratory linear models with different term
#' specifications and recommends the simplest model that adequately captures
#' data structure.
#'
#' Uses adjusted R-squared, information criteria (AIC/BIC), and residual
#' patterns to balance model fit against complexity. Recommendations consider:
#' \itemize{
#'   \item Linear effects: whether factors have substantial linear relationships
#'   \item Interactions: whether factor interactions improve fit significantly
#'   \item Quadratic effects: whether curvature is present and non-negligible
#' }
#'
#' @param data Data frame containing response and factor columns.
#' @param response Name of the response variable.
#' @param factor_names Character vector of factor names.
#' @param criterion Information criterion to use for model comparison.
#'   One of \code{"AIC"} (default), \code{"BIC"}, or \code{"adjR2"}
#'   (adjusted R-squared).
#' @param alpha Significance threshold for model comparison. Default is 0.05.
#' @param verbose Logical; if \code{TRUE}, prints detailed analysis of each
#'   model term. Default is \code{FALSE}.
#'
#' @return A list with class \code{brsm_term_recommendation} containing:
#'   \itemize{
#'     \item \code{recommended_terms}: Recommended \code{model_terms}
#'       specification (character string).
#'     \item \code{fit_comparison}: Data frame with one row per model term,
#'       showing number of parameters, R-squared, adjusted R-squared, AIC,
#'       BIC, and residual standard error.
#'     \item \code{diagnostics}: List of diagnostic results including
#'       linear_rsq (R-squared from linear-only model), has_interaction_value
#'       (whether interactions improve fit substantially), has_quadratic_value
#'       (whether quadratics improve fit substantially), mean_abs_residual
#'       (mean absolute residual from recommended model).
#'     \item \code{call}: The function call.
#'   }
#'   The result has a custom \code{print} method that displays recommendations
#'   and diagnostics clearly.
#'
#' @examples
#' \dontrun{
#' # Linear data: will recommend first_order
#' dat_linear <- data.frame(
#'   x1 = runif(50, -1, 1),
#'   x2 = runif(50, -1, 1),
#'   y = 5 + 2*x1 + 3*x2 + rnorm(50, sd = 0.5)
#' )
#' recommend_brsm_model_terms(dat_linear, "y", c("x1", "x2"))
#'
#' # Quadratic data: will recommend second_order
#' dat_quad <- data.frame(
#'   x1 = runif(100, -1, 1),
#'   x2 = runif(100, -1, 1),
#'   y = 5 + 2*x1 - x1^2 + 0.5*x1*x2 + rnorm(100, sd = 0.3)
#' )
#' recommend_brsm_model_terms(dat_quad, "y", c("x1", "x2"))
#' }
#'
#' @export
recommend_brsm_model_terms <- function(data,
                                      response,
                                      factor_names,
                                      criterion = c("AIC", "BIC", "adjR2"),
                                      alpha = 0.05,
                                      verbose = FALSE) {
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
  criterion <- match.arg(criterion)

  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha < 0 || alpha > 1) {
    stop("alpha must be a number between 0 and 1.")
  }
  if (!is.logical(verbose) || length(verbose) != 1L) {
    stop("verbose must be TRUE or FALSE.")
  }

  # Build and fit models with different term specifications
  model_list <- list(
    first_order = .brsm_build_exploratory_model(
      data, response, factor_names, "first_order"
    ),
    first_order_twi = .brsm_build_exploratory_model(
      data, response, factor_names, "first_order_twi"
    ),
    pure_quadratic = .brsm_build_exploratory_model(
      data, response, factor_names, "pure_quadratic"
    ),
    second_order = .brsm_build_exploratory_model(
      data, response, factor_names, "second_order"
    )
  )

  # Extract fit statistics
  fit_comparison <- .brsm_extract_fit_statistics(model_list)

  # Run diagnostics
  diagnostics <- .brsm_compute_term_diagnostics(
    data = data,
    response = response,
    factor_names = factor_names,
    model_list = model_list,
    fit_comparison = fit_comparison,
    alpha = alpha,
    verbose = verbose
  )

  # Determine recommendation
  recommended <- .brsm_determine_recommendation(
    fit_comparison = fit_comparison,
    diagnostics = diagnostics,
    criterion = criterion,
    verbose = verbose
  )

  # Print summary
  cat("\nModel Term Recommendation:\n\n")
  cat("Recommended model_terms: ", recommended, "\n\n")
  print(fit_comparison, row.names = FALSE)
  cat("\n")

  if (diagnostics$has_interaction_value && recommended %in% c("second_order", "first_order_twi")) {
    cat("✓ Two-way interactions improve fit significantly\n")
  }
  if (diagnostics$has_quadratic_value && recommended %in% c("second_order", "pure_quadratic")) {
    cat("✓ Quadratic terms improve fit significantly\n")
  }
  if (!diagnostics$has_interaction_value && !diagnostics$has_quadratic_value) {
    cat("✓ Linear model is adequate; no evidence for higher-order terms\n")
  }
  cat("\n")

  result <- list(
    recommended_terms = recommended,
    fit_comparison = fit_comparison,
    diagnostics = diagnostics,
    call = match.call()
  )

  class(result) <- c("brsm_term_recommendation", "list")
  invisible(result)
}


#' @export
print.brsm_term_recommendation <- function(x, ...) {
  cat("\nModel Term Recommendation:\n\n")
  cat("Recommended: ", x$recommended_terms, "\n\n")

  cat("Fit Comparison:\n")
  print(x$fit_comparison, row.names = FALSE)
  cat("\n")

  invisible(x)
}


# Build an exploratory lm model with specified term structure
.brsm_build_exploratory_model <- function(data, response, factor_names, model_terms) {
  linear_terms <- factor_names

  interaction_terms <- if (model_terms %in% c("second_order", "first_order_twi") &&
                            length(factor_names) > 1) {
    pairs <- utils::combn(factor_names, 2, simplify = FALSE)
    vapply(pairs, function(p) paste0(p[1], ":", p[2]), character(1))
  } else character(0)

  quadratic_terms <- if (model_terms %in% c("second_order", "pure_quadratic")) {
    paste0("I(", factor_names, "^2)")
  } else character(0)

  all_terms <- c(linear_terms, interaction_terms, quadratic_terms)
  formula_text <- paste(response, "~", paste(all_terms, collapse = " + "))
  model_formula <- stats::as.formula(formula_text)

  stats::lm(model_formula, data = data)
}


# Extract and summarize fit statistics from models
.brsm_extract_fit_statistics <- function(model_list) {
  stat_names <- c("model_terms", "n_params", "n_obs", "r_squared", "adj_r_squared",
                  "AIC", "BIC", "sigma")

  stats_list <- lapply(names(model_list), function(model_name) {
    fit <- model_list[[model_name]]
    if (inherits(fit, "try-error")) {
      return(data.frame(
        model_terms = model_name,
        n_params = NA_real_,
        n_obs = NA_real_,
        r_squared = NA_real_,
        adj_r_squared = NA_real_,
        AIC = NA_real_,
        BIC = NA_real_,
        sigma = NA_real_
      ))
    }

    summary_fit <- summary(fit)
    n_params <- length(stats::coef(fit))
    n_obs <- nrow(fit$model)

    data.frame(
      model_terms = model_name,
      n_params = n_params,
      n_obs = n_obs,
      r_squared = summary_fit$r.squared,
      adj_r_squared = summary_fit$adj.r.squared,
      AIC = stats::AIC(fit),
      BIC = stats::BIC(fit),
      sigma = summary_fit$sigma,
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, stats_list)
}


# Compute diagnostics for term value (interactions, quadratics, etc.)
.brsm_compute_term_diagnostics <- function(data, response, factor_names, model_list,
                                            fit_comparison, alpha, verbose) {
  fit_linear <- model_list$first_order
  fit_interaction <- model_list$first_order_twi
  fit_quad <- model_list$pure_quadratic
  fit_second <- model_list$second_order

  # Check for interaction value: does first_order_twi improve over first_order?
  has_interaction <- FALSE
  if (!inherits(fit_interaction, "try-error") && !inherits(fit_linear, "try-error")) {
    anova_result <- stats::anova(fit_linear, fit_interaction)
    p_value_interaction <- anova_result[2, "Pr(>F)"]
    has_interaction <- !is.na(p_value_interaction) && p_value_interaction < alpha

    if (verbose && !is.na(p_value_interaction)) {
      cat("Interaction test p-value:", signif(p_value_interaction, 3), "\n")
    }
  }

  # Check for quadratic value: does pure_quadratic improve over first_order?
  has_quadratic <- FALSE
  if (!inherits(fit_quad, "try-error") && !inherits(fit_linear, "try-error")) {
    anova_result <- stats::anova(fit_linear, fit_quad)
    p_value_quadratic <- anova_result[2, "Pr(>F)"]
    has_quadratic <- !is.na(p_value_quadratic) && p_value_quadratic < alpha

    if (verbose && !is.na(p_value_quadratic)) {
      cat("Quadratic test p-value:", signif(p_value_quadratic, 3), "\n")
    }
  }

  # Residual analysis on recommended model
  linear_rsq <- summary(fit_linear)$r.squared
  recommended_fit <- if (has_quadratic || has_interaction) {
    fit_second
  } else {
    fit_linear
  }

  mean_abs_residual <- mean(abs(stats::residuals(recommended_fit)), na.rm = TRUE)

  list(
    linear_rsq = linear_rsq,
    has_interaction_value = has_interaction,
    has_quadratic_value = has_quadratic,
    mean_abs_residual = mean_abs_residual,
    n_factors = length(factor_names)
  )
}


# Determine which model to recommend based on fit statistics
.brsm_determine_recommendation <- function(fit_comparison, diagnostics, criterion, verbose) {
  # Remove any NA rows
  fit_clean <- fit_comparison[!is.na(fit_comparison$AIC), ]

  if (nrow(fit_clean) == 0) {
    warning("All model fits failed; defaulting to second_order.")
    return("second_order")
  }

  # Score models based on criterion
  if (criterion == "AIC") {
    best_idx <- which.min(fit_clean$AIC)
  } else if (criterion == "BIC") {
    best_idx <- which.min(fit_clean$BIC)
  } else if (criterion == "adjR2") {
    best_idx <- which.max(fit_clean$adj_r_squared)
  }

  best_model <- fit_clean$model_terms[best_idx]

  # Apply parsimony check: if within 2 AIC points of best and simpler,
  # prefer simpler model
  best_aic <- fit_clean$AIC[best_idx]
  complexity_order <- c("first_order", "first_order_twi", "pure_quadratic", "second_order")

  for (model_name in complexity_order) {
    model_idx <- which(fit_clean$model_terms == model_name)
    if (length(model_idx) == 0) next

    model_aic <- fit_clean$AIC[model_idx]
    if (model_aic <= best_aic + 2) {
      best_model <- model_name
      break
    }
  }

  if (verbose) {
    cat("\nModel selection based on", criterion, ":\n")
    print(fit_clean[, c("model_terms", criterion, "n_params")])
  }

  best_model
}