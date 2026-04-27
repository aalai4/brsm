#' Identify Stationary Points of Response Surface
#'
#' Computes stationary (critical) points of the response surface from posterior
#' draws, including computation of the eigenvalues of the Hessian matrix for
#' classification.
#'
#' @param object A `brsm_fit` object, `brmsfit` object, or
#'   data frame of posterior draws with Bayesian coefficient columns.
#' @param factor_names Character vector of factor names
#'   (if object is data frame).
#' @param kappa_thresh Threshold for computing condition numbers of Hessian.
#' @param diagnostics Diagnostic reporting level. One of \code{"none"}
#'   (default), \code{"basic"}, or \code{"full"}. When not \code{"none"},
#'   a diagnostics list is attached as an attribute on the returned data frame.
#' @param auto_guidance Logical; if \code{TRUE} (default), emit guidance
#'   warnings when exclusion rates are high or threshold sensitivity is unstable.
#' @param sensitivity_thresholds Numeric vector of condition-number thresholds
#'   used when \code{diagnostics = "full"}. Defaults to
#'   \code{c(1e8, 1e10, 1e12)}.
#'
#' @return A data frame with stationary points and their properties.
#'
#' @export
stationary_point <- function(object,
                             factor_names = NULL,
                             kappa_thresh = 1e10,
                             diagnostics = c("none", "basic", "full"),
                             auto_guidance = TRUE,
                             sensitivity_thresholds = c(1e8, 1e10, 1e12)) {
  UseMethod("stationary_point")
}

#' @rdname stationary_point
#' @export
stationary_point.brsm_fit <- function(
    object,
    factor_names = NULL,
    kappa_thresh = 1e10,
    diagnostics = c("none", "basic", "full"),
    auto_guidance = TRUE,
    sensitivity_thresholds = c(1e8, 1e10, 1e12)) {
  if (is.null(factor_names)) {
    factor_names <- object$factor_names
  }
  stationary_point.default(
    as_brsm_draws(object),
    factor_names = factor_names,
    kappa_thresh = kappa_thresh,
    diagnostics = diagnostics,
    auto_guidance = auto_guidance,
    sensitivity_thresholds = sensitivity_thresholds
  )
}

#' @rdname stationary_point
#' @export
stationary_point.default <- function(
    object,
    factor_names = NULL,
    kappa_thresh = 1e10,
    diagnostics = c("none", "basic", "full"),
    auto_guidance = TRUE,
    sensitivity_thresholds = c(1e8, 1e10, 1e12)) {
  draws <- .brsm_validate_draws(object)
  if (is.null(factor_names)) {
    stop("factor_names must be supplied if object does not contain metadata.")
  }
  factor_names <- .brsm_validate_factor_names(factor_names)

  diagnostics <- match.arg(diagnostics)
  if (!is.numeric(kappa_thresh) || length(kappa_thresh) != 1L ||
      !is.finite(kappa_thresh) || kappa_thresh <= 0) {
    stop("kappa_thresh must be a finite positive numeric scalar.")
  }
  if (!is.logical(auto_guidance) || length(auto_guidance) != 1L ||
      is.na(auto_guidance)) {
    stop("auto_guidance must be TRUE or FALSE.")
  }
  if (!is.numeric(sensitivity_thresholds) || length(sensitivity_thresholds) == 0L ||
      any(!is.finite(sensitivity_thresholds)) || any(sensitivity_thresholds <= 0)) {
    stop("sensitivity_thresholds must be a non-empty finite positive numeric vector.")
  }
  sensitivity_thresholds <- sort(unique(as.numeric(sensitivity_thresholds)))

  n_draws <- nrow(draws)
  n_factors <- length(factor_names)

  # Store column names once
  draw_cols <- colnames(draws)

  # Linear coefficients
  linear_cols <- paste0("b_", factor_names)

  if (!all(linear_cols %in% draw_cols)) {
    stop(
      "missing linear columns: ",
      paste(setdiff(linear_cols, draw_cols), collapse = ", ")
    )
  }

  if (!all(vapply(
    draws[, linear_cols, drop = FALSE],
    is.numeric,
    logical(1)
  ))) {
    stop("draws must contain numeric columns for linear terms.")
  }

  b <- as.matrix(draws[, linear_cols, drop = FALSE])
  storage.mode(b) <- "numeric"

  # Hessian matrices
  H <- .brsm_hessian_array(draws, factor_names)

  solve_detail <- stationary_points_batch_details(
    h_array = H,
    b_matrix = b,
    kappa_thresh = as.numeric(kappa_thresh)
  )
  x_star <- solve_detail$x_star
  colnames(x_star) <- factor_names

  status_code <- as.integer(solve_detail$status_code)
  kappa_proxy <- as.numeric(solve_detail$kappa_proxy)
  status_levels <- .brsm_stationary_status_levels()
  status_label <- .brsm_stationary_status_labels(status_code)

  excluded <- rowSums(is.na(x_star)) > 0
  n_excluded <- sum(excluded)
  pct_excluded <- n_excluded / n_draws

  if (n_excluded > 0) {
    warning(n_excluded, " draws have near-singular Hessians and were set to NA.")
  }

  x_star <- as.data.frame(x_star)
  rownames(x_star) <- NULL

  if (diagnostics != "none") {
    diag_info <- list(
      n_draws = n_draws,
      n_factors = n_factors,
      kappa_thresh = as.numeric(kappa_thresh),
      n_excluded = n_excluded,
      pct_excluded = pct_excluded,
      status_code = status_code,
      status_label = status_label,
      kappa_proxy = kappa_proxy,
      status_counts = .brsm_stationary_status_counts(status_code),
      function_name = "stationary_point"
    )

    if (diagnostics == "full") {
      thresholds <- sort(unique(c(sensitivity_thresholds, as.numeric(kappa_thresh))))
      excl_counts <- integer(length(thresholds))
      status_by_threshold <- vector("list", length(thresholds))

      for (i in seq_along(thresholds)) {
        detail_i <- stationary_points_batch_details(
          h_array = H,
          b_matrix = b,
          kappa_thresh = thresholds[i]
        )
        status_i <- as.integer(detail_i$status_code)
        counts_i <- .brsm_stationary_status_counts(status_i)
        excl_counts[i] <- sum(status_i != 0L)
        status_by_threshold[[i]] <- data.frame(
          threshold = thresholds[i],
          ok = unname(counts_i[["ok"]]),
          lapack_fail = unname(counts_i[["lapack_fail"]]),
          invalid_lu_diag = unname(counts_i[["invalid_lu_diag"]]),
          kappa_exceeded = unname(counts_i[["kappa_exceeded"]]),
          n_excluded = excl_counts[i],
          pct_excluded = excl_counts[i] / n_draws,
          stringsAsFactors = FALSE
        )
      }

      sensitivity_df <- data.frame(
        threshold = thresholds,
        n_excluded = excl_counts,
        pct_excluded = excl_counts / n_draws,
        stringsAsFactors = FALSE
      )
      excl_range_pp <- 100 * (max(sensitivity_df$pct_excluded) -
        min(sensitivity_df$pct_excluded))
      stability_flag <- if (excl_range_pp <= 5) {
        "stable"
      } else if (excl_range_pp <= 15) {
        "moderately_sensitive"
      } else {
        "unstable"
      }

      diag_info$exclusion_by_threshold <- sensitivity_df
      diag_info$status_by_threshold <- do.call(rbind, status_by_threshold)
      diag_info$stability_flag <- stability_flag
      diag_info$exclusion_range_pp <- excl_range_pp
    }

    attr(x_star, "diagnostics") <- diag_info

    if (isTRUE(auto_guidance)) {
      if (pct_excluded > 0.2) {
        warning(
          sprintf(
            "High stationary-point exclusion rate (%.1f%%). ",
            100 * pct_excluded
          ),
          "Consider checking sensitivity at kappa_thresh = 1e12 and using ridge analysis."
        )
      } else if (pct_excluded > 0.15) {
        warning(
          sprintf(
            "Moderately high stationary-point exclusion rate (%.1f%%). ",
            100 * pct_excluded
          ),
          "Consider a sensitivity check across kappa thresholds."
        )
      }

      if (diagnostics == "full" &&
          identical(attr(x_star, "diagnostics")$stability_flag, "unstable")) {
        warning(
          "Stationary-point exclusions are unstable across sensitivity thresholds; ",
          "treat this as a low-information geometry warning."
        )
      }
    }
  }

  x_star
}