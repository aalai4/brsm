#' Check MCMC Diagnostics for a brsm Fit
#'
#' Provides a compact post-fit diagnostics summary focused on fixed-effect
#' posterior parameters: Rhat, bulk ESS, tail ESS, and NUTS diagnostics
#' (divergences, treedepth saturation, BFMI).
#'
#' @param object A \code{brsm_fit} object from [fit_brsm()] or a
#'   \code{brmsfit} object.
#' @param rhat_threshold Threshold above which Rhat is flagged.
#' @param ess_bulk_min Minimum recommended bulk ESS.
#' @param ess_tail_min Minimum recommended tail ESS.
#' @param treedepth_limit Optional treedepth saturation threshold. If
#'   \code{NULL}, attempts to infer from fit control settings and falls back
#'   to \code{10}.
#' @param bfmi_threshold BFMI threshold below which chains are flagged.
#' @param verbose Logical; if \code{TRUE}, prints a concise summary.
#'
#' @return A list with components:
#'   \code{overview} (one-row data frame), \code{parameters} (per-parameter
#'   diagnostics), and \code{passed} (logical scalar).
#'
#' @examples
#' \dontrun{
#' fit <- fit_brsm(dat, response = "y", factor_names = c("x1", "x2"))
#' diag <- check_brsm_fit(fit)
#' diag$overview
#' }
#' @export
.brsm_extract_fit <- function(object, caller = "function") {
  fit <- object
  if (inherits(object, "brsm_fit")) {
    fit <- object$fit
  }

  if (!inherits(fit, "brmsfit")) {
    stop("object must be a brsm_fit or brmsfit object in ", caller, "().")
  }

  fit
}

.brsm_get_max_treedepth <- function(fit) {
  td <- NA_integer_

  td <- tryCatch(
    {
      as.integer(fit$fit@sim$control$max_treedepth)
    },
    error = function(e) NA_integer_
  )

  if (!is.na(td)) {
    return(td)
  }

  td <- tryCatch(
    {
      if (is.function(fit$fit$metadata)) {
        md <- fit$fit$metadata()
        if ("max_treedepth" %in% names(md)) {
          return(as.integer(md$max_treedepth))
        }
      }
      NA_integer_
    },
    error = function(e) NA_integer_
  )

  td
}

.brsm_compute_bfmi <- function(np) {
  if (is.null(np) || !is.data.frame(np) ||
    !all(c("Chain", "Parameter", "Value") %in% names(np))) {
    return(data.frame(chain = integer(0), bfmi = numeric(0)))
  }

  e_rows <- np[np$Parameter == "energy__", c("Chain", "Value"), drop = FALSE]
  if (nrow(e_rows) == 0) {
    return(data.frame(chain = integer(0), bfmi = numeric(0)))
  }

  split_e <- split(e_rows$Value, e_rows$Chain)
  bfmi <- vapply(split_e, function(en) {
    if (length(en) < 2 || stats::var(en) <= 0) {
      return(NA_real_)
    }
    mean(diff(en)^2) / stats::var(en)
  }, numeric(1))

  data.frame(
    chain = as.integer(names(bfmi)),
    bfmi = as.numeric(bfmi),
    stringsAsFactors = FALSE
  )
}

check_brsm_fit <- function(object,
                           rhat_threshold = 1.01,
                           ess_bulk_min = 400,
                           ess_tail_min = 400,
                           treedepth_limit = NULL,
                           bfmi_threshold = 0.3,
                           verbose = TRUE) {
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("package 'brms' is required for check_brsm_fit().")
  }

  fit <- .brsm_extract_fit(object, caller = "check_brsm_fit")

  if (!is.null(treedepth_limit) &&
    (!is.numeric(treedepth_limit) ||
      length(treedepth_limit) != 1 ||
      !is.finite(treedepth_limit))) {
    stop("treedepth_limit must be NULL or a finite numeric scalar.")
  }
  if (!is.numeric(bfmi_threshold) ||
    length(bfmi_threshold) != 1 ||
    !is.finite(bfmi_threshold)) {
    stop("bfmi_threshold must be a finite numeric scalar.")
  }

  summ <- brms::summary(fit)
  fixed <- summ$fixed

  if (is.null(fixed) || nrow(fixed) == 0) {
    stop("No fixed-effect diagnostics found in fit summary.")
  }

  fixed <- as.data.frame(fixed)
  fixed$parameter <- rownames(fixed)
  rownames(fixed) <- NULL

  rhat_col <- if ("Rhat" %in% names(fixed)) "Rhat" else NULL
  ess_bulk_col <- if ("Bulk_ESS" %in% names(fixed)) {
    "Bulk_ESS"
  } else if ("Eff.Sample" %in% names(fixed)) {
    "Eff.Sample"
  } else {
    NULL
  }
  ess_tail_col <- if ("Tail_ESS" %in% names(fixed)) "Tail_ESS" else NULL

  diag_df <- data.frame(
    parameter = fixed$parameter,
    rhat = if (!is.null(rhat_col)) as.numeric(fixed[[rhat_col]]) else NA_real_,
    ess_bulk = if (!is.null(ess_bulk_col)) {
      as.numeric(fixed[[ess_bulk_col]])
    } else {
      NA_real_
    },
    ess_tail = if (!is.null(ess_tail_col)) {
      as.numeric(fixed[[ess_tail_col]])
    } else {
      NA_real_
    },
    stringsAsFactors = FALSE
  )

  n_rhat_bad <- if (all(is.na(diag_df$rhat))) {
    NA_integer_
  } else {
    sum(diag_df$rhat > rhat_threshold, na.rm = TRUE)
  }
  n_ess_bulk_bad <- if (all(is.na(diag_df$ess_bulk))) {
    NA_integer_
  } else {
    sum(diag_df$ess_bulk < ess_bulk_min, na.rm = TRUE)
  }
  n_ess_tail_bad <- if (all(is.na(diag_df$ess_tail))) {
    NA_integer_
  } else {
    sum(diag_df$ess_tail < ess_tail_min, na.rm = TRUE)
  }

  np <- tryCatch(brms::nuts_params(fit), error = function(e) NULL)
  divergences <- NA_integer_
  divergence_rate <- NA_real_
  n_max_treedepth_hits <- NA_integer_
  treedepth_hit_rate <- NA_real_
  bfmi_min <- NA_real_
  n_bfmi_below_threshold <- NA_integer_
  if (!is.null(np) && is.data.frame(np) &&
    all(c("Parameter", "Value") %in% names(np))) {
    div_rows <- np[np$Parameter == "divergent__", , drop = FALSE]
    if (nrow(div_rows) > 0) {
      divergences <- sum(div_rows$Value > 0, na.rm = TRUE)
      divergence_rate <- divergences / nrow(div_rows)
    }

    if (is.null(treedepth_limit)) {
      treedepth_limit <- .brsm_get_max_treedepth(fit)
      if (is.na(treedepth_limit)) {
        treedepth_limit <- 10
      }
    }

    td_rows <- np[np$Parameter == "treedepth__", , drop = FALSE]
    if (nrow(td_rows) > 0) {
      n_max_treedepth_hits <- sum(
        td_rows$Value >= treedepth_limit,
        na.rm = TRUE
      )
      treedepth_hit_rate <- n_max_treedepth_hits / nrow(td_rows)
    }

    bfmi_df <- .brsm_compute_bfmi(np)
    if (nrow(bfmi_df) > 0) {
      bfmi_min <- min(bfmi_df$bfmi, na.rm = TRUE)
      n_bfmi_below_threshold <- sum(bfmi_df$bfmi < bfmi_threshold, na.rm = TRUE)
    }
  }

  rhat_ok <- is.na(n_rhat_bad) || n_rhat_bad == 0
  ess_bulk_ok <- is.na(n_ess_bulk_bad) || n_ess_bulk_bad == 0
  ess_tail_ok <- is.na(n_ess_tail_bad) || n_ess_tail_bad == 0
  divergence_ok <- is.na(divergences) || divergences == 0
  treedepth_ok <- is.na(n_max_treedepth_hits) || n_max_treedepth_hits == 0
  bfmi_ok <- is.na(n_bfmi_below_threshold) || n_bfmi_below_threshold == 0
  passed <- rhat_ok && ess_bulk_ok && ess_tail_ok &&
    divergence_ok && treedepth_ok && bfmi_ok

  overview <- data.frame(
    n_parameters = nrow(diag_df),
    rhat_max = if (all(is.na(diag_df$rhat))) {
      NA_real_
    } else {
      max(diag_df$rhat, na.rm = TRUE)
    },
    n_rhat_over_threshold = n_rhat_bad,
    ess_bulk_min_observed = if (all(is.na(diag_df$ess_bulk))) {
      NA_real_
    } else {
      min(diag_df$ess_bulk, na.rm = TRUE)
    },
    n_ess_bulk_below_min = n_ess_bulk_bad,
    ess_tail_min_observed = if (all(is.na(diag_df$ess_tail))) {
      NA_real_
    } else {
      min(diag_df$ess_tail, na.rm = TRUE)
    },
    n_ess_tail_below_min = n_ess_tail_bad,
    divergences = divergences,
    divergence_rate = divergence_rate,
    treedepth_limit = treedepth_limit,
    n_max_treedepth_hits = n_max_treedepth_hits,
    treedepth_hit_rate = treedepth_hit_rate,
    bfmi_threshold = bfmi_threshold,
    bfmi_min = bfmi_min,
    n_bfmi_below_threshold = n_bfmi_below_threshold,
    passed = passed,
    stringsAsFactors = FALSE
  )

  if (isTRUE(verbose)) {
    message(
      "check_brsm_fit: passed=", passed,
      ", rhat_max=", format(overview$rhat_max, digits = 4),
      ", min_bulk_ess=", format(overview$ess_bulk_min_observed, digits = 4),
      ", min_tail_ess=", format(overview$ess_tail_min_observed, digits = 4),
      ", divergences=", overview$divergences,
      ", treedepth_hits=", overview$n_max_treedepth_hits,
      ", bfmi_min=", format(overview$bfmi_min, digits = 4)
    )
  }

  list(
    overview = overview,
    parameters = diag_df,
    passed = passed
  )
}
