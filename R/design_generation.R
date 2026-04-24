#' Generate Standard Response Surface Designs
#'
#' Generates coded experimental designs commonly used in response surface
#' methodology (RSM): Central Composite Design (CCD) and Box-Behnken Design
#' (BBD).
#'
#' Designs are returned in coded units by default, where factor levels are
#' centered at 0 and typically span approximately \eqn{[-1, 1]}.
#'
#' @param factor_names Character vector of factor names.
#' @param design Design type. One of \code{"ccd"} or \code{"bbd"}.
#' @param n_center Number of center-point replicates.
#' @param alpha Axial distance used for CCD. Either \code{"rotatable"},
#'   \code{"orthogonal"}, or a positive numeric scalar.
#' @param randomize Logical; if \code{TRUE}, randomize run order.
#' @param seed Optional random seed used when \code{randomize = TRUE}.
#' @param output Output scale. One of \code{"coded"} (default) or
#'   \code{"natural"}.
#' @param ranges Optional named list of factor ranges in natural units, each a
#'   numeric vector of length 2. Required when \code{output = "natural"}.
#'
#' @return A data frame of design runs with columns:
#'   \itemize{
#'   \item factor columns
#'   \item \code{.design_type} (ccd/bbd)
#'   \item \code{.block} (factorial/axial/center for CCD; edge/center for BBD)
#'   \item \code{.run_id} (sequential run index)
#'   }
#'   If \code{output = "coded"} and \code{ranges} are supplied, coding
#'   metadata is attached as \code{brsm_coding} for downstream decode helpers.
#'
#' @export
generate_brsm_design <- function(factor_names,
                                 design = c("ccd", "bbd"),
                                 n_center = 4,
                                 alpha = c("rotatable", "orthogonal"),
                                 randomize = TRUE,
                                 seed = NULL,
                                 output = c("coded", "natural"),
                                 ranges = NULL) {
  factor_names <- .brsm_validate_factor_names(factor_names)
  design <- match.arg(design)
  output <- match.arg(output)

  if (!is.numeric(n_center) || length(n_center) != 1L || !is.finite(n_center) ||
      n_center < 0) {
    stop("n_center must be a finite numeric scalar >= 0.")
  }
  n_center <- as.integer(n_center)

  if (!is.logical(randomize) || length(randomize) != 1L || is.na(randomize)) {
    stop("randomize must be TRUE or FALSE.")
  }

  k <- length(factor_names)

  if (design == "bbd" && k < 3L) {
    stop("BBD requires at least 3 factors.")
  }

  if (output == "natural" && is.null(ranges)) {
    stop("ranges must be provided when output = 'natural'.")
  }

  if (design == "ccd") {
    if (is.character(alpha)) {
      if (length(alpha) > 1L) {
        alpha <- alpha[[1L]]
      }
      alpha_key <- tolower(alpha)
      if (!alpha_key %in% c("rotatable", "orthogonal")) {
        stop("alpha must be 'rotatable', 'orthogonal', or a positive numeric scalar.")
      }
      alpha <- alpha_key
    }

    coded <- .brsm_generate_ccd(
      factor_names = factor_names,
      n_center = n_center,
      alpha = alpha
    )
  } else {
    coded <- .brsm_generate_bbd(
      factor_names = factor_names,
      n_center = n_center
    )
  }

  if (isTRUE(randomize) && nrow(coded) > 1L) {
    if (!is.null(seed)) {
      set.seed(seed)
    }
    coded <- coded[sample.int(nrow(coded)), , drop = FALSE]
  }

  coded$.run_id <- seq_len(nrow(coded))

  if (output == "coded") {
    if (!is.null(ranges)) {
      coding <- .brsm_ranges_to_coding(ranges, factor_names)
      attr(coded, "brsm_coding") <- coding
      class(coded) <- unique(c("brsm_coded_data", class(coded)))
    }
    return(coded)
  }

  natural <- .brsm_decode_from_ranges(coded, ranges, factor_names)
  natural$.design_type <- coded$.design_type
  natural$.block <- coded$.block
  natural$.run_id <- coded$.run_id
  natural
}


.brsm_generate_ccd <- function(factor_names, n_center, alpha) {
  k <- length(factor_names)

  alpha_val <- .brsm_resolve_alpha(alpha, k)

  # Full 2^k factorial at -1/+1
  fac_grid <- as.matrix(expand.grid(rep(list(c(-1, 1)), k), KEEP.OUT.ATTRS = FALSE))
  colnames(fac_grid) <- factor_names
  fac_df <- as.data.frame(fac_grid)
  fac_df$.design_type <- "ccd"
  fac_df$.block <- "factorial"

  # Axial points at +/- alpha on each axis
  axial <- matrix(0, nrow = 2 * k, ncol = k)
  for (i in seq_len(k)) {
    axial[(2 * i) - 1, i] <- -alpha_val
    axial[2 * i, i] <- alpha_val
  }
  colnames(axial) <- factor_names
  axial_df <- as.data.frame(axial)
  axial_df$.design_type <- "ccd"
  axial_df$.block <- "axial"

  # Center points
  center_df <- as.data.frame(matrix(0, nrow = n_center, ncol = k))
  colnames(center_df) <- factor_names
  center_df$.design_type <- "ccd"
  center_df$.block <- "center"

  out <- rbind(fac_df, axial_df)
  if (n_center > 0L) {
    out <- rbind(out, center_df)
  }
  rownames(out) <- NULL
  out
}


.brsm_generate_bbd <- function(factor_names, n_center) {
  k <- length(factor_names)
  pairs <- utils::combn(seq_len(k), 2, simplify = FALSE)

  edge_rows <- list()
  for (pair in pairs) {
    signs <- as.matrix(expand.grid(c(-1, 1), c(-1, 1), KEEP.OUT.ATTRS = FALSE))
    for (r in seq_len(nrow(signs))) {
      x <- rep(0, k)
      x[pair[1]] <- signs[r, 1]
      x[pair[2]] <- signs[r, 2]
      edge_rows[[length(edge_rows) + 1L]] <- x
    }
  }

  edge_mat <- do.call(rbind, edge_rows)
  colnames(edge_mat) <- factor_names
  edge_df <- as.data.frame(edge_mat)
  edge_df$.design_type <- "bbd"
  edge_df$.block <- "edge"

  center_df <- as.data.frame(matrix(0, nrow = n_center, ncol = k))
  colnames(center_df) <- factor_names
  center_df$.design_type <- "bbd"
  center_df$.block <- "center"

  out <- edge_df
  if (n_center > 0L) {
    out <- rbind(out, center_df)
  }
  rownames(out) <- NULL
  out
}


.brsm_resolve_alpha <- function(alpha, k) {
  if (is.numeric(alpha) && length(alpha) == 1L && is.finite(alpha) && alpha > 0) {
    return(as.numeric(alpha))
  }

  if (is.character(alpha) && length(alpha) == 1L) {
    key <- tolower(alpha)
    if (key == "rotatable") {
      return((2^k)^(1 / 4))
    }
    if (key == "orthogonal") {
      return(sqrt(k))
    }
  }

  stop("alpha must be 'rotatable', 'orthogonal', or a positive numeric scalar.")
}


.brsm_ranges_to_coding <- function(ranges, factor_names) {
  if (!is.list(ranges) || is.null(names(ranges))) {
    stop("ranges must be a named list.")
  }
  if (!all(factor_names %in% names(ranges))) {
    stop("ranges must include all factor_names.")
  }

  factor_coding <- vector("list", length(factor_names))
  names(factor_coding) <- factor_names

  for (f in factor_names) {
    r <- ranges[[f]]
    if (!is.numeric(r) || length(r) != 2L || any(!is.finite(r))) {
      stop("Each range must be a finite numeric vector of length 2.")
    }
    r <- sort(as.numeric(r))
    center <- mean(r)
    scale <- diff(r) / 2
    if (!is.finite(scale) || scale <= 0) {
      stop("Range for factor '", f, "' must have positive width.")
    }

    factor_coding[[f]] <- list(
      center = center,
      scale = scale,
      method = "range",
      formula = paste0(
        f, "_coded = (", f, " - ", signif(center, 8), ") / ", signif(scale, 8)
      ),
      decode_formula = paste0(
        f, " = ", f, "_coded * ", signif(scale, 8), " + ", signif(center, 8)
      ),
      original_range = r,
      coded_range = c(-1, 1)
    )
  }

  list(method = "range", factors = factor_coding)
}


.brsm_decode_from_ranges <- function(coded_df, ranges, factor_names) {
  coding <- .brsm_ranges_to_coding(ranges, factor_names)

  out <- coded_df
  for (f in factor_names) {
    center <- coding$factors[[f]]$center
    scale <- coding$factors[[f]]$scale
    out[[f]] <- as.numeric(out[[f]]) * scale + center
  }

  out
}