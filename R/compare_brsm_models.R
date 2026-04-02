#' Compare Bayesian Models for brsm
#'
#' Computes model comparison criteria across multiple fitted Bayesian models.
#'
#' @param models A named list of \code{brsm_fit} and/or \code{brmsfit}
#'   objects.
#' @param criterion Comparison criterion: \code{"loo"} or \code{"waic"}.
#' @param ... Additional arguments passed to \code{brms::loo()} or
#'   \code{brms::waic()}.
#'
#' @return A list with components:
#'   \code{criterion}, \code{estimates} (named list of criterion objects), and
#'   \code{comparison} (data frame from \code{loo_compare}).
#'
#' @examples
#' \dontrun{
#' fit1 <- fit_brsm(dat, response = "y", factor_names = c("x1", "x2"))
#' fit2 <- fit_brsm(dat,
#'   response = "y", factor_names = c("x1", "x2"),
#'   prior = brms::prior(normal(0, 1), class = "b")
#' )
#' cmp <- compare_brsm_models(
#'   list(default = fit1, stronger_prior = fit2), "loo"
#' )
#' cmp$comparison
#' }
#' @export
compare_brsm_models <- function(models,
                                criterion = c("loo", "waic"),
                                ...) {
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("package 'brms' is required for compare_brsm_models().")
  }

  criterion <- match.arg(criterion)

  if (!is.list(models) || length(models) < 2) {
    stop("models must be a list with at least two brsm_fit/brmsfit objects.")
  }

  if (is.null(names(models)) || any(names(models) == "")) {
    names(models) <- paste0("model_", seq_along(models))
  }

  fits <- lapply(
    models,
    function(m) .brsm_extract_fit(m, caller = "compare_brsm_models")
  )

  estimates <- if (criterion == "loo") {
    lapply(fits, function(f) brms::loo(f, ...))
  } else {
    lapply(fits, function(f) brms::waic(f, ...))
  }

  cmp <- do.call(brms::loo_compare, unname(estimates))
  cmp_df <- as.data.frame(cmp)
  cmp_df$model <- rownames(cmp)
  rownames(cmp_df) <- NULL

  cmp_df <- cmp_df[, c("model", setdiff(names(cmp_df), "model")), drop = FALSE]

  list(
    criterion = criterion,
    estimates = estimates,
    comparison = cmp_df
  )
}