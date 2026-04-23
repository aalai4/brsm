#' Prepare Factor Columns for BRSM Modeling
#'
#' Applies lightweight centering/scaling to factor columns for response-surface
#' modeling and stores coding metadata needed to decode predictions later.
#'
#' @param data Data frame containing factor columns.
#' @param factor_names Character vector of factor column names to code.
#' @param method Coding method. One of \code{"zscore"} (mean/sd) or
#'   \code{"range"} (midpoint/half-range; approximately maps observed range to
#'   \code{[-1, 1]}), or \code{"identity"} (no transform; records explicit
#'   identity coding metadata for already-coded predictors).
#'
#' @return A data frame with coded factor columns. The returned object has class
#'   \code{brsm_coded_data} and contains attribute \code{brsm_coding}, a list
#'   with formulas and scaling parameters for each factor.
#'
#' @examples
#' dat <- data.frame(
#'   x1 = rnorm(10, 5, 2), x2 = runif(10, 10, 20), y = rnorm(10)
#' )
#' coded <- prepare_brsm_data(
#'   dat, factor_names = c("x1", "x2"), method = "range"
#' )
#' str(attr(coded, "brsm_coding"))
#' @export
prepare_brsm_data <- function(data,
                              factor_names,
                              method = c("zscore", "range", "identity")) {
  if (!is.data.frame(data)) {
    stop("data must be a data.frame.")
  }

  factor_names <- .brsm_validate_factor_names(factor_names)
  method <- match.arg(method)

  missing_factors <- setdiff(factor_names, names(data))
  if (length(missing_factors) > 0) {
    stop(
      "Missing factor columns in data: ",
      paste(missing_factors, collapse = ", ")
    )
  }

  if (!all(vapply(data[factor_names], is.numeric, logical(1)))) {
    stop("All factor columns must be numeric.")
  }

  out <- data
  factor_coding <- vector("list", length(factor_names))
  names(factor_coding) <- factor_names

  for (f in factor_names) {
    x <- as.numeric(out[[f]])

    if (all(is.na(x))) {
      stop("Factor '", f, "' contains only NA values.")
    }

    if (method == "zscore") {
      center <- mean(x, na.rm = TRUE)
      scale <- stats::sd(x, na.rm = TRUE)
    } else if (method == "range") {
      x_range <- range(x, na.rm = TRUE)
      center <- mean(x_range)
      scale <- diff(x_range) / 2
    } else {
      center <- 0
      scale <- 1
    }

    if (method != "identity" && (!is.finite(scale) || scale <= 0)) {
      stop(
        "Factor '", f, "' has zero/invalid scale under method '", method,
        "'. Provide a factor with variation."
      )
    }

    out[[f]] <- (x - center) / scale

    factor_coding[[f]] <- list(
      center = center,
      scale = scale,
      method = method,
      formula = paste0(
        f, "_coded = (", f, " - ",
        signif(center, 8), ") / ", signif(scale, 8)
      ),
      decode_formula = paste0(
        f, " = ", f, "_coded * ",
        signif(scale, 8), " + ", signif(center, 8)
      ),
      original_range = range(x, na.rm = TRUE),
      coded_range = range(out[[f]], na.rm = TRUE)
    )
  }

  coding <- list(
    method = method,
    factors = factor_coding
  )

  attr(out, "brsm_coding") <- coding
  class(out) <- unique(c("brsm_coded_data", class(out)))
  out
}


#' Extract Stored brsm Coding Metadata
#'
#' Retrieves coding metadata from a \code{brsm_coded_data} object or a
#' \code{brsm_fit} object created from coded data.
#'
#' @param object A \code{brsm_coded_data} data frame or \code{brsm_fit} object.
#'
#' @return A coding list with \code{method} and \code{factors} entries.
#' @export
get_brsm_coding <- function(object) {
  if (inherits(object, "brsm_fit")) {
    coding <- object$coding
  } else {
    coding <- attr(object, "brsm_coding", exact = TRUE)
  }

  if (is.null(coding)) {
    stop("No brsm coding metadata found on object.")
  }

  if (!is.list(coding) || is.null(coding$factors) || !is.list(coding$factors)) {
    stop("Invalid brsm coding metadata format.")
  }

  coding
}


#' Decode Coded Factor Columns Back to Original Scale
#'
#' Uses stored coding metadata to decode factor columns back to original units.
#' This is useful for decoding prediction grids or summaries generated on coded
#' scales.
#'
#' @param data Data frame with coded factor columns.
#' @param coding Optional coding metadata from [get_brsm_coding()]. If
#'   \code{NULL}, the function uses \code{attr(data, "brsm_coding")},
#'   if present.
#' @param factor_names Optional character vector of factors to decode. If
#'   \code{NULL}, decodes all factors available in \code{coding}.
#' @param overwrite Logical; if \code{TRUE}, decoded values overwrite existing
#'   factor columns. If \code{FALSE}, decoded columns are appended with
#'   \code{suffix}.
#' @param suffix Suffix appended to decoded column names when
#'   \code{overwrite = FALSE}.
#'
#' @return A data frame with decoded factor values.
#'
#' @examples
#' dat <- data.frame(
#'   x1 = rnorm(10, 5, 2), x2 = runif(10, 10, 20), y = rnorm(10)
#' )
#' coded <- prepare_brsm_data(dat, c("x1", "x2"), method = "range")
#' decoded <- decode_brsm_data(coded)
#' @export
decode_brsm_data <- function(data,
                             coding = NULL,
                             factor_names = NULL,
                             overwrite = TRUE,
                             suffix = "_decoded") {
  if (!is.data.frame(data)) {
    stop("data must be a data.frame.")
  }

  if (is.null(coding)) {
    coding <- attr(data, "brsm_coding", exact = TRUE)
  }
  if (is.null(coding)) {
    stop("coding is NULL and no 'brsm_coding' attribute found on data.")
  }

  if (!is.list(coding) || is.null(coding$factors) || !is.list(coding$factors)) {
    stop("coding must be a list with a 'factors' list entry.")
  }

  available_factors <- names(coding$factors)
  if (is.null(factor_names)) {
    factor_names <- available_factors
  }
  factor_names <- .brsm_validate_factor_names(factor_names)

  missing_meta <- setdiff(factor_names, available_factors)
  if (length(missing_meta) > 0) {
    stop(
      "No coding metadata found for factors: ",
      paste(missing_meta, collapse = ", ")
    )
  }

  out <- data
  for (f in factor_names) {
    if (!f %in% names(out)) {
      stop("data is missing coded factor column: ", f)
    }

    center <- coding$factors[[f]]$center
    scale <- coding$factors[[f]]$scale

    decoded <- as.numeric(out[[f]]) * scale + center

    if (overwrite) {
      out[[f]] <- decoded
    } else {
      out[[paste0(f, suffix)]] <- decoded
    }
  }

  if (overwrite) {
    class(out) <- setdiff(class(out), "brsm_coded_data")
    attr(out, "brsm_coding") <- NULL
  }

  out
}