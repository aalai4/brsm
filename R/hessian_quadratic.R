hessian_quadratic <- function(draws, factor_names) {
  draws <- .brsm_validate_draws(draws)
  factor_names <- .brsm_validate_factor_names(factor_names)

  hess_array <- .brsm_hessian_array(draws, factor_names)
  hess_df <- expand.grid(
    draw = seq_len(dim(hess_array)[1]),
    row_factor = factor_names,
    col_factor = factor_names,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  hess_df$value <- as.vector(hess_array)

  hess_df
}
