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

  if (isTRUE(include_intercept)) {
    prior_list[[length(prior_list) + 1L]] <- brms::set_prior(
      .brsm_prior_string(coefficient_family, 0, intercept_scale, student_df),
      class = "Intercept"
    )
  }

  need_global_b <- !is.null(available_b_coefs) && any(lengths(resolved$unmatched) > 0L)
  if (need_global_b) {
    prior_list[[length(prior_list) + 1L]] <- brms::set_prior(
      .brsm_prior_string(coefficient_family, 0, linear_scale, student_df),
      class = "b"
    )
  }

  for (f in resolved$matched$linear) {
    prior_list[[length(prior_list) + 1L]] <- brms::set_prior(
      .brsm_prior_string(coefficient_family, 0, linear_scale, student_df),
      class = "b",
      coef = f
    )
  }

  if (length(resolved$matched$interaction) > 0L) {
    for (coef_name in resolved$matched$interaction) {
      prior_list[[length(prior_list) + 1L]] <- brms::set_prior(
        .brsm_prior_string(coefficient_family, 0, interaction_scale, student_df),
        class = "b",
        coef = coef_name
      )
    }
  }

  if (length(resolved$matched$quadratic) > 0L) {
    for (coef_name in resolved$matched$quadratic) {
      prior_list[[length(prior_list) + 1L]] <- brms::set_prior(
        .brsm_prior_string(coefficient_family, 0, quadratic_scale, student_df),
        class = "b",
        coef = coef_name
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


.brsm_default_prior_b_coefs <- function(response,
                                        linear_terms,
                                        interaction_terms,
                                        quadratic_terms,
                                        include_interactions,
                                        include_quadratic,
                                        data) {
  if (is.null(data) || !is.data.frame(data) ||
      is.null(response) || !is.character(response) || length(response) != 1L ||
      !response %in% names(data)) {
    return(NULL)
  }

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

  formula_text <- paste(response, "~", paste(rhs_terms, collapse = " + "))
  model_formula <- stats::as.formula(formula_text)

  prior_df <- tryCatch(
    as.data.frame(brms::default_prior(model_formula, data = data)),
    error = function(e) NULL
  )

  if (is.null(prior_df) || !all(c("class", "coef") %in% names(prior_df))) {
    return(NULL)
  }

  coefs <- unique(prior_df$coef[prior_df$class == "b"])
  coefs <- coefs[!is.na(coefs) & nzchar(coefs)]
  if (length(coefs) == 0L) NULL else as.character(coefs)
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