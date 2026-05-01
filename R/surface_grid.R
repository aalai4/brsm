#' Generate a Regular Grid over Factor Ranges
#'
#' Creates an expanded prediction grid spanning the specified factor ranges,
#' which is used internally for surface evaluation and contour plotting.
#'
#' @param ranges A named list of numeric vectors of length 2, one per factor,
#'   giving the \code{c(min, max)} range for each factor.
#' @param n Integer; number of equally-spaced grid points per factor.
#'   Must be \eqn{\geq 2}. Default is \code{50}.
#' @param center Optional named numeric vector. If supplied, this point is
#'   guaranteed to appear in the grid (appended if not already present).
#'
#' @return A data frame with one column per factor and \code{n^p} rows
#'   (where \code{p} is the number of factors), representing the full
#'   factorial grid of evaluation points.
#'
#' @export
surface_grid <- function(ranges,
                         n = 50,
                         center = NULL) {
  if (is.null(ranges)) {
    stop("ranges must be supplied.")
  }

  if (!is.list(ranges)) {
    stop("ranges must be a named list.")
  }

  if (is.null(names(ranges))) {
    stop("ranges must be a named list with factor names.")
  }

  if (!is.numeric(n) || length(n) != 1 || n < 2 || !is.finite(n)) {
    stop("n must be a finite integer >= 2.")
  }

  n <- as.integer(n)
  factor_names <- names(ranges)

  grid_vals <- vector("list", length(factor_names))
  names(grid_vals) <- factor_names

  for (f in factor_names) {
    r <- ranges[[f]]

    if (!is.numeric(r) || length(r) != 2) {
      stop(paste0("range for factor '", f, "' must be numeric of length 2."))
    }

    if (any(!is.finite(r))) {
      stop(paste0("range for factor '", f, "' must contain finite values."))
    }

    r <- sort(r)

    grid_vals[[f]] <- seq(r[1], r[2], length.out = n)
  }

  grid <- expand.grid(grid_vals, KEEP.OUT.ATTRS = FALSE)

  if (!is.null(center)) {
    if (!is.numeric(center)) {
      stop("center must be a numeric vector.")
    }

    if (is.null(names(center))) {
      stop("center must be a named numeric vector.")
    }

    if (!all(factor_names %in% names(center))) {
      stop("center must contain values for all factors.")
    }

    if (length(center) != length(factor_names)) {
      stop("center must have the same length as the number of factors.")
    }

    center_row <- as.data.frame(t(center[factor_names]))

    is_dup <- any(vapply(
      seq_len(nrow(grid)),
      function(i) {
        row <- as.numeric(grid[i, factor_names])
        all(abs(row - center[factor_names]) < 1e-8)
      },
      logical(1)
    ))

    if (!is_dup) {
      grid <- rbind(center_row, grid)
    }
  }

  rownames(grid) <- NULL

  grid
}
