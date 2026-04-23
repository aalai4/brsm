.brsm_validate_factor_names <- function(factor_names,
                                        require_two = FALSE,
                                        two_factors_message = NULL) {
  if (is.null(factor_names)) {
    stop("factor_names must be supplied.")
  }

  factor_names <- as.character(factor_names)

  if (length(factor_names) == 0) {
    stop("factor_names must contain at least one factor.")
  }

  if (require_two && length(factor_names) != 2) {
    if (!is.null(two_factors_message)) {
      stop(two_factors_message)
    }
    stop("factor_names must contain exactly two factors.")
  }

  factor_names
}

.brsm_validate_draws <- function(
  draws,
  empty_message = "draws must contain at least one posterior draw."
) {
  draws <- as.data.frame(draws)

  if (nrow(draws) == 0) {
    stop(empty_message)
  }

  draws
}

.brsm_validate_probs <- function(
  probs,
  require_length = NULL,
  require_finite = FALSE,
  message = NULL
) {
  bad <- !is.numeric(probs)

  if (!bad && !is.null(require_length)) {
    bad <- length(probs) != require_length
  }

  if (!bad && require_finite) {
    bad <- any(!is.finite(probs))
  }

  if (!bad) {
    bad <- any(probs < 0) || any(probs > 1)
  }

  if (bad) {
    if (!is.null(message)) {
      stop(message)
    }
    if (require_finite) {
      stop("probs must contain finite values between 0 and 1.")
    }
    stop("probs must contain values between 0 and 1.")
  }

  sort(probs)
}

.brsm_validate_bayesian_input <- function(object) {
  # Ensure input is Bayesian-only
  if (inherits(object, "brsm_fit")) {
    if (is.null(object$fit) || !inherits(object$fit, "brmsfit")) {
      stop("brsm_fit object must contain a valid brmsfit model in `$fit`.")
    }
    return(invisible(object))
  }

  if (inherits(object, "brmsfit")) {
    # brmsfit is acceptable
    return(invisible(object))
  }

  if (is.data.frame(object)) {
    # Data frame is acceptable if it has posterior draw columns (b_*)
    b_cols <- grep("^b_", names(object), value = TRUE)
    if (length(b_cols) == 0) {
      stop(
        "Input data frame lacks Bayesian posterior columns. ",
        "Data frames must have columns named 'b_Intercept', 'b_x1', etc. ",
        "Consider: (1) using a brmsfit object from brms::brm(), or ",
        "(2) converting your posterior draws to the required format."
      )
    }
    return(invisible(object))
  }

  # Reject non-Bayesian input types
  obj_class <- paste(class(object), collapse = ", ")
  stop(
    "brsm workflow input must be strictly Bayesian. ",
    "Accepted input types: (1) brsm_fit object from fit_brsm(), ",
    "(2) brmsfit object from brms::brm(), or ",
    "(3) data frame with Bayesian posterior coefficient columns ",
    "(b_Intercept, b_x1, etc.). ",
    "Received object of class: ", obj_class, ". ",
    "This package does not support frequentist models."
  )
}

.brsm_enforce_coding_policy <- function(data,
                                        factor_names,
                                        coding,
                                        coding_policy = c(
                                          "warn", "error", "ignore"
                                        ),
                                        caller = "fit_brsm") {
  coding_policy <- match.arg(coding_policy)

  if (!is.null(coding)) {
    return(invisible(coding_policy))
  }

  means <- vapply(
    data[factor_names],
    function(x) mean(x, na.rm = TRUE),
    numeric(1)
  )
  sds <- vapply(
    data[factor_names],
    function(x) stats::sd(x, na.rm = TRUE),
    numeric(1)
  )

  centered <- all(abs(means) < 0.1)
  scaled <- all(abs(sds - 1) < 0.1)

  msg <- paste0(
    caller, "() called with no brsm coding metadata. ",
    "This means prepare_brsm_data() was not used (or metadata was lost). ",
    "Uncoded predictors can cause numerically unstable quadratic ",
    "Bayesian fits. ",
    "Recommended: run prepare_brsm_data(data, factor_names, ",
    "method='zscore', 'range', or 'identity') before fitting."
  )

  if (centered && scaled) {
    msg <- paste0(
      msg,
      " Predictors appear approximately centered/scaled already; ",
      "if this was intentional, set coding_policy='ignore' to ",
      "silence this message."
    )
  }

  if (coding_policy == "error") {
    stop(msg)
  }

  if (coding_policy == "warn") {
    warning(msg)
  }

  invisible(coding_policy)
}

.brsm_check_columns <- function(required, obj, message) {
  if (!all(required %in% names(obj))) {
    stop(message)
  }
}

# Resolve an interaction column by trying both factor orderings and both
# separator styles (`:` from brms; `.` from data.frame check.names coercion).
# Returns the first matching column name, or NULL if none is found.
.brsm_find_interaction_col <- function(f1, f2, draw_cols) {
  candidates <- c(
    paste0("b_", f1, ":", f2),
    paste0("b_", f2, ":", f1),
    paste0("b_", f1, ".", f2),
    paste0("b_", f2, ".", f1)
  )
  found <- intersect(candidates, draw_cols)
  if (length(found) == 0L) NULL else found[[1L]]
}

# Resolve a quadratic column across canonical and sanitized naming
# variants produced by brms/posterior/data.frame coercion.
.brsm_find_quadratic_col <- function(f, draw_cols) {
  canonical <- paste0("b_I(", f, "^2)")
  candidates <- unique(c(
    canonical,
    paste0("b_I", f, "E2"),
    paste0("b_I", f, "E2."),
    paste0("b_I.", f, ".2."),
    paste0("b_I.", f, ".2"),
    paste0("b_I", f, ".2."),
    paste0("b_I", f, ".2"),
    make.names(canonical)
  ))

  found <- intersect(candidates, draw_cols)
  if (length(found) > 0L) {
    return(found[[1L]])
  }

  f_escaped <- gsub("([][{}()+*^$|\\\\?.])", "\\\\\\\\\\1", f)
  idx <- grep(
    paste0("^b_?I.*", f_escaped, ".*(\\\\^2|E2|\\\\.2\\\\.?|2\\\\)?)$"),
    draw_cols,
    perl = TRUE
  )

  if (length(idx) == 0L) NA_character_ else draw_cols[idx[[1L]]]
}

.brsm_hessian_array <- function(draws, factor_names) {
  draws <- .brsm_validate_draws(draws)
  factor_names <- .brsm_validate_factor_names(factor_names)

  quad_cols <- vapply(
    factor_names,
    function(f) .brsm_find_quadratic_col(f, colnames(draws)),
    character(1)
  )
  missing_idx <- which(is.na(quad_cols) | quad_cols == "")
  missing_cols <- paste0("b_I(", factor_names[missing_idx], "^2)")

  n_factors <- length(factor_names)
  n_draws <- nrow(draws)

  if (length(missing_idx) > 0) {
    stop("missing quadratic columns: ", paste(missing_cols, collapse = ", "))
  }
  if (!all(vapply(draws[, quad_cols, drop = FALSE], is.numeric, logical(1)))) {
    stop("draws must contain numeric columns for quadratic terms.")
  }
  quad <- as.matrix(draws[, quad_cols, drop = FALSE])

  hess_array <- array(0, dim = c(n_draws, n_factors, n_factors))

  for (i in seq_len(n_factors)) {
    hess_array[, i, i] <- 2 * quad[, i]
  }

  if (n_factors > 1) {
    for (i in seq_len(n_factors - 1)) {
      for (j in (i + 1):n_factors) {
        inter_col <- .brsm_find_interaction_col(
          factor_names[i], factor_names[j], colnames(draws)
        )
        coef <- if (!is.null(inter_col)) draws[[inter_col]] else rep(0, n_draws)
        if (!is.numeric(coef)) {
          stop("draws must contain numeric columns for interaction terms.")
        }
        hess_array[, i, j] <- coef
        hess_array[, j, i] <- coef
      }
    }
  }

  dimnames(hess_array) <- list(
    draw = NULL,
    row = factor_names,
    col = factor_names
  )

  hess_array
}