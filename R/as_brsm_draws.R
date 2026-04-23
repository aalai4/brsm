#' Convert Model Output to brsm Draws Format
#'
#' Converts posterior draws to the brsm-compatible format with
#' standardized column names (b_Intercept, b_x1, b_I(x1^2), etc.).
#'
#' @param object A `brsm_fit` object, `brmsfit` object, or a data frame of
#'   posterior draws/coefficients.
#' @param factor_names Character vector of factor names (e.g., c("x1", "x2")).
#' @param ... Additional arguments (currently unused).
#'
#' @return A data frame with brsm-compatible column names:
#'   - b_Intercept: intercept
#'   - b_<factor>: linear term for each factor
#'   - b_I(<factor>^2): quadratic term for each factor
#'   - b_<factor>:<factor>: interaction terms
#'
#' @examples
#' # From posterior draws in a data frame
#' draws_raw <- data.frame(
#'   b_Intercept = rnorm(100, 5, 0.2),
#'   b_x1 = rnorm(100, 2, 0.1),
#'   b_x2 = rnorm(100, 4, 0.1),
#'   "b_I(x1^2)" = rnorm(100, -1, 0.05),
#'   "b_I(x2^2)" = rnorm(100, -2, 0.05),
#'   "b_x1:x2" = rnorm(100, 0.5, 0.05),
#'   check.names = FALSE
#' )
#' draws <- as_brsm_draws(draws_raw, factor_names = c("x1", "x2"))
#' head(draws)
#'
#' @export
as_brsm_draws <- function(object, factor_names, ...) {
  UseMethod("as_brsm_draws")
}

#' @rdname as_brsm_draws
#' @export
as_brsm_draws.brsm_fit <- function(object, factor_names = NULL, ...) {
  if (is.null(object$fit) || !inherits(object$fit, "brmsfit")) {
    stop("brsm_fit object must contain a valid brmsfit model in `$fit`.")
  }
  if (is.null(factor_names)) {
    factor_names <- object$factor_names
  }
  as_brsm_draws.brmsfit(object$fit, factor_names = factor_names, ...)
}

#' @rdname as_brsm_draws
#' @export
as_brsm_draws.brmsfit <- function(object, factor_names = NULL, ...) {
  if (is.null(factor_names)) {
    stop("factor_names must be supplied for brmsfit objects.")
  }
  draws_raw <- as.data.frame(object, ...)
  draw_cols <- names(draws_raw)

  required_main <- c("b_Intercept", paste0("b_", factor_names))
  missing_cols <- setdiff(required_main, draw_cols)

  # Quadratic term names must follow canonical brsm naming: b_I(x1^2).
  quadratic_missing <- character(0)
  for (f in factor_names) {
    quad_col <- .brsm_find_quadratic_col(f, draw_cols)
    if (is.na(quad_col) || quad_col == "") {
      quadratic_missing <- c(quadratic_missing, paste0("b_I(", f, "^2)"))
    }
  }

  missing_interactions <- character(0)
  if (length(factor_names) > 1) {
    for (i in seq_len(length(factor_names) - 1)) {
      for (j in (i + 1):length(factor_names)) {
        f1 <- factor_names[i]
        f2 <- factor_names[j]
        interaction_candidates <- c(
          paste0("b_", f1, ":", f2),
          paste0("b_", f2, ":", f1),
          paste0("b_", f1, ".", f2),
          paste0("b_", f2, ".", f1)
        )
        if (!any(interaction_candidates %in% draw_cols)) {
          missing_interactions <- c(
            missing_interactions,
            paste0("b_", f1, ":", f2)
          )
        }
      }
    }
  }

  any_missing <- length(missing_cols) > 0 ||
    length(quadratic_missing) > 0 ||
    length(missing_interactions) > 0
  if (any_missing) {
    stop(
      "brms fit is missing required posterior coefficient columns: ",
      paste(
        c(missing_cols, quadratic_missing, missing_interactions),
        collapse = ", "
      )
    )
  }

  draws <- as_brsm_draws.data.frame(draws_raw, factor_names = factor_names)
  rownames(draws) <- NULL
  draws
}

#' @rdname as_brsm_draws
#' @export
as_brsm_draws.data.frame <- function(object, factor_names, ...) {
  # Assume data frame is already in the right format or needs column renaming
  # Try to detect and rename if necessary
  df <- object

  # Normalize accepted quadratic naming variants to canonical b_I(x^2).
  df_renamed <- df
  for (f in factor_names) {
    canonical_col <- paste0("b_I(", f, "^2)")
    quad_col <- .brsm_find_quadratic_col(f, names(df_renamed))
    if (!is.na(quad_col) && quad_col != "" &&
      quad_col != canonical_col &&
      !(canonical_col %in% names(df_renamed))) {
      names(df_renamed)[names(df_renamed) == quad_col] <- canonical_col
    }
  }

  # Check if already in brsm format
  required_cols <- c(
    "b_Intercept",
    paste0("b_", factor_names),
    paste0("b_I(", factor_names, "^2)")
  )

  # Detect and include interaction columns among factor_names.
  interaction_cols <- character(0)
  if (length(factor_names) > 1) {
    for (i in seq_len(length(factor_names) - 1)) {
      for (j in (i + 1):length(factor_names)) {
        f1 <- factor_names[i]
        f2 <- factor_names[j]
        candidate_cols <- c(
          paste0("b_", f1, ":", f2),
          paste0("b_", f2, ":", f1),
          paste0("b_", f1, ".", f2),
          paste0("b_", f2, ".", f1)
        )
        found <- intersect(candidate_cols, names(df_renamed))
        if (length(found) > 0) {
          canonical_interaction <- paste0("b_", f1, ":", f2)
          interaction_col <- found[[1]]
          if (interaction_col != canonical_interaction &&
            !(canonical_interaction %in% names(df_renamed))) {
            names(df_renamed)[
              names(df_renamed) == interaction_col
            ] <- canonical_interaction
            interaction_col <- canonical_interaction
          }
          interaction_cols <- c(interaction_cols, interaction_col)
        }
      }
    }
  }

  # Keep additional coefficient columns that already follow b_* naming.
  # This preserves higher-order or custom terms without manual mapping.
  extra_b_cols <- grep("^b_", names(df_renamed), value = TRUE)
  keep_cols <- unique(c(required_cols, interaction_cols, extra_b_cols))

  # Return only brsm-format columns
  brsm_cols <- intersect(keep_cols, names(df_renamed))
  missing_required <- setdiff(required_cols, brsm_cols)

  if (length(missing_required) > 0) {
    stop(
      "Data frame does not contain recognizable brsm column names. ",
      "Missing required columns: ",
      paste(missing_required, collapse = ", ")
    )
  }

  df_renamed[, brsm_cols, drop = FALSE]
}

#' @rdname as_brsm_draws
#' @export
as_brsm_draws.default <- function(object, factor_names, ...) {
  stop(
    "as_brsm_draws() does not support objects of class: ",
    paste(class(object), collapse = ", "),
    "\nProvide a brsm_fit, brmsfit, or data frame with ",
    "brsm-compatible coefficient columns."
  )
}