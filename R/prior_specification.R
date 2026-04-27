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
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("package 'brms' is required for specify_brsm_priors().")
  }

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

  prior_list <- list()

  if (isTRUE(include_intercept)) {
    prior_list[[length(prior_list) + 1L]] <- brms::set_prior(
      .brsm_prior_string(coefficient_family, 0, intercept_scale, student_df),
      class = "Intercept"
    )
  }

  for (f in factor_names) {
    prior_list[[length(prior_list) + 1L]] <- brms::set_prior(
      .brsm_prior_string(coefficient_family, 0, linear_scale, student_df),
      class = "b",
      coef = f
    )
  }

  include_interactions <- model_terms %in% c("second_order", "first_order_twi")
  if (include_interactions && length(factor_names) > 1L) {
    pairs <- utils::combn(factor_names, 2, simplify = FALSE)
    for (pair in pairs) {
      prior_list[[length(prior_list) + 1L]] <- brms::set_prior(
        .brsm_prior_string(coefficient_family, 0, interaction_scale, student_df),
        class = "b",
        coef = paste0(pair[1], ":", pair[2])
      )
    }
  }

  include_quadratic <- model_terms %in% c("second_order", "pure_quadratic")
  if (include_quadratic) {
    for (f in factor_names) {
      prior_list[[length(prior_list) + 1L]] <- brms::set_prior(
        .brsm_prior_string(coefficient_family, 0, quadratic_scale, student_df),
        class = "b",
        coef = paste0("I(", f, "^2)")
      )
    }
  }

  if (isTRUE(include_sigma)) {
    prior_list[[length(prior_list) + 1L]] <- brms::set_prior(
      paste0("student_t(3, 0, ", signif(sigma_prior_scale, 6), ")"),
      class = "sigma"
    )
  }

  do.call(c, prior_list)
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