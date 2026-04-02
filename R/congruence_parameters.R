#' Compute Congruence-Specific Transformed Parameters
#'
#' Computes the standard Response Surface Analysis (RSA) congruence hypothesis
#' parameters (a1 through a5) from posterior draws. These transformations map
#' coefficients from a two-factor quadratic model to domain-specific RSA
#' constructs.
#'
#' The congruence parameters follow the Edwards & Parry (1993) / Shanock et al.
#' (2010) RSA framework:
#'   - \code{a1}: Sum of linear coefficients (b1 + b2)
#'   - \code{a2}: Sum of curvature coefficients (b11 + b22)
#'   - \code{a3}: Interaction term (b12)
#'   - \code{a4}: Curvature interaction (-2 * b11 * b22)
#'   - \code{a5}: Curvature difference (|b11 - b22|)
#'
#' @param object A \code{brsm_fit} object from [fit_brsm()], a \code{brmsfit}
#'   object from \code{brms::brm()}, or a data frame of posterior draws with
#'   required Bayesian coefficient columns.
#' @param factor_names Character vector of exactly two factor names.
#'   If \code{NULL}, will attempt to infer from object metadata.
#'
#' @return A data frame with 5 columns (\code{a1} through \code{a5}), one row
#'   per posterior draw. Each draw is a complete realization of the RSA
#'   congruence hypothesis parameters.
#'
#' @examples
#' \dontrun{
#' # After fitting a model:
#' fit <- fit_brsm(data, response = "y", factor_names = c("x1", "x2"))
#' cong_params <- congruence_parameters(fit, factor_names = c("x1", "x2"))
#' head(cong_params)
#' }
#'
#' @export
congruence_parameters <- function(object, factor_names = NULL) {
  UseMethod("congruence_parameters")
}

#' @rdname congruence_parameters
#' @export
congruence_parameters.brsm_fit <- function(object, factor_names = NULL) {
  if (is.null(factor_names)) {
    factor_names <- object$factor_names
  }
  congruence_parameters.default(
    as_brsm_draws(object),
    factor_names = factor_names
  )
}

#' @rdname congruence_parameters
#' @export
congruence_parameters.brmsfit <- function(object, factor_names = NULL) {
  congruence_parameters.default(
    as.data.frame(object),
    factor_names = factor_names
  )
}

#' @rdname congruence_parameters
#' @export
congruence_parameters.default <- function(object, factor_names = NULL) {
  # Extract Bayesian input
  draws <- .brsm_validate_draws(object)

  if (is.null(factor_names)) {
    stop("factor_names must be supplied if object does not contain metadata.")
  }

  factor_names <- .brsm_validate_factor_names(
    factor_names,
    require_two = TRUE,
    two_factors_message = "congruence_parameters requires exactly two factors."
  )

  f1 <- factor_names[1]
  f2 <- factor_names[2]

  # Build expected column names for two-factor quadratic model
  b1_col <- paste0("b_", f1)
  b2_col <- paste0("b_", f2)
  b11_col <- paste0("b_I(", f1, "^2)")
  b22_col <- paste0("b_I(", f2, "^2)")
  b12_col <- paste0("b_", f1, ":", f2)

  b11_col_use <- .brsm_find_quadratic_col(f1, colnames(draws))
  if (is.na(b11_col_use)) {
    stop(
      "Could not find quadratic term column for ", f1,
      ". Expected: ", b11_col
    )
  }

  b22_col_use <- .brsm_find_quadratic_col(f2, colnames(draws))
  if (is.na(b22_col_use)) {
    stop(
      "Could not find quadratic term column for ", f2,
      ". Expected: ", b22_col
    )
  }

  if (b12_col %in% colnames(draws)) {
    b12_col_use <- b12_col
  } else {
    stop(
      "Could not find interaction term column for ", f1, ":", f2,
      ". Tried: ", b12_col
    )
  }

  # Validate required columns exist
  required_cols <- c(b1_col, b2_col, b11_col_use, b22_col_use, b12_col_use)
  missing_cols <- setdiff(required_cols, colnames(draws))
  if (length(missing_cols) > 0) {
    stop(
      "Missing coefficient columns: ", paste(missing_cols, collapse = ", "),
      ". Expected for quadratic model with factors: ",
      paste(factor_names, collapse = ", ")
    )
  }

  # Extract coefficient vectors
  b1 <- draws[[b1_col]]
  b2 <- draws[[b2_col]]
  b11 <- draws[[b11_col_use]]
  b22 <- draws[[b22_col_use]]
  b12 <- draws[[b12_col_use]]

  # Compute RSA congruence parameters
  a1 <- b1 + b2 # Sum of linear coefficients
  a2 <- b11 + b22 # Sum of curvature coefficients
  a3 <- b12 # Interaction term
  a4 <- -2 * b11 * b22 # Curvature interaction (surface interaction)
  a5 <- abs(b11 - b22) # Curvature difference

  data.frame(
    a1 = a1,
    a2 = a2,
    a3 = a3,
    a4 = a4,
    a5 = a5
  )
}


#' Summarize Congruence Parameters with Credible Intervals
#'
#' Computes posterior summaries (mean, median, credible intervals) of the
#' five RSA congruence hypothesis parameters (a1 through a5).
#'
#' @param object A \code{brsm_fit} object from [fit_brsm()], a \code{brmsfit}
#'   object from \code{brms::brm()}, or a data frame of posterior draws with
#'   required Bayesian coefficient columns.
#' @param factor_names Character vector of exactly two factor names.
#'   If \code{NULL}, will attempt to infer from object metadata.
#' @param probs Numeric vector of credible interval probabilities. Default is
#'   \code{c(0.025, 0.25, 0.5, 0.75, 0.975)} for 95% and 50% intervals.
#'
#' @return A data frame with columns:
#'   \code{parameter} (a1-a5), \code{mean}, \code{median}, and one column per
#'   probability in \code{probs} (named \code{q_XX} where XX is the percentile).
#'
#' @examples
#' \dontrun{
#' fit <- fit_brsm(data, response = "y", factor_names = c("x1", "x2"))
#' summary_cong <- summarize_congruence(fit, factor_names = c("x1", "x2"))
#' print(summary_cong)
#' }
#'
#' @export
summarize_congruence <- function(object,
                                 factor_names = NULL,
                                 probs = c(0.025, 0.25, 0.5, 0.75, 0.975)) {
  # Get congruence parameters
  cong_params <- congruence_parameters(object, factor_names = factor_names)

  probs <- .brsm_validate_probs(probs)

  # Compute summaries for each parameter
  param_names <- colnames(cong_params)
  summaries <- lapply(param_names, function(pname) {
    x <- cong_params[[pname]]
    mean_val <- mean(x, na.rm = TRUE)
    median_val <- stats::median(x, na.rm = TRUE)
    quantiles <- stats::quantile(x, probs = probs, na.rm = TRUE)

    # Format quantile column names
    q_names <- paste0("q_", sprintf("%03d", as.integer(probs * 100)))

    data.frame(
      parameter = pname,
      mean = mean_val,
      median = median_val,
      as.list(quantiles),
      check.names = FALSE,
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  })

  # Bind summaries
  result <- do.call(rbind, summaries)
  rownames(result) <- NULL

  # Rename generic quantile columns to match specified probs
  q_names_expected <- paste0("q_", sprintf("%03d", as.integer(probs * 100)))
  new_col_names <- c("parameter", "mean", "median", q_names_expected)
  colnames(result) <- new_col_names[1:ncol(result)]

  result
}


#' Apply ROPE (Region of Practical Equivalence) to Congruence Parameters
#'
#' Tests whether RSA congruence hypothesis parameters fall within a specified
#' Region of Practical Equivalence (ROPE) and computes posterior proportions
#' inside and outside the ROPE.
#'
#' The ROPE is application-dependent. Common values in RSA studies are
#' \code{c(-0.1, 0.1)} or domain-specific bounds based on
#' effect size conventions.
#'
#' @param object A \code{brsm_fit} object from [fit_brsm()], a \code{brmsfit}
#'   object from \code{brms::brm()}, or a data frame of posterior draws with
#'   required Bayesian coefficient columns.
#' @param factor_names Character vector of exactly two factor names.
#'   If \code{NULL}, will attempt to infer from object metadata.
#' @param rope Numeric vector of length 2 specifying the lower and upper bounds
#'   of the ROPE. Default is \code{c(-0.1, 0.1)}.
#'
#' @return A data frame with columns:
#'   \code{parameter} (a1-a5), \code{inside_rope}, \code{outside_rope},
#'   and \code{decision} (a string classifying the evidence).
#'
#'   Decision rules (approximate):
#'   - "Practical Equivalence": >97.5\% of posterior is inside ROPE
#'   - "Not Equivalent": >97.5\% of posterior is outside ROPE
#'   - "Inconclusive": Otherwise
#'
#' @examples
#' \dontrun{
#' fit <- fit_brsm(data, response = "y", factor_names = c("x1", "x2"))
#' rope_results <- rope_congruence(
#'   fit,
#'   factor_names = c("x1", "x2"),
#'   rope = c(-0.1, 0.1)
#' )
#' print(rope_results)
#' }
#'
#' @export
rope_congruence <- function(object,
                            factor_names = NULL,
                            rope = c(-0.1, 0.1)) {
  # Validate ROPE specification
  if (!is.numeric(rope) || length(rope) != 2) {
    stop("rope must be a numeric vector of length 2 (lower, upper bounds).")
  }
  if (rope[1] >= rope[2]) {
    stop("rope[1] must be less than rope[2].")
  }

  # Get congruence parameters
  cong_params <- congruence_parameters(object, factor_names = factor_names)

  param_names <- colnames(cong_params)
  rope_results <- lapply(param_names, function(pname) {
    x <- cong_params[[pname]]

    # Classify each draw relative to ROPE
    inside <- (x >= rope[1] & x <= rope[2])
    outside <- (x < rope[1] | x > rope[2])

    # Compute proportions
    prob_inside <- mean(inside, na.rm = TRUE)
    prob_outside <- mean(outside, na.rm = TRUE)

    # Decision rule (adapted from Bayesian equivalence testing)
    if (prob_inside > 0.975) {
      decision <- "Practical Equivalence"
    } else if (prob_outside > 0.975) {
      decision <- "Not Equivalent"
    } else {
      decision <- "Inconclusive"
    }

    data.frame(
      parameter = pname,
      inside_rope = prob_inside,
      outside_rope = prob_outside,
      decision = decision,
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  })

  result <- do.call(rbind, rope_results)
  rownames(result) <- NULL
  result
}