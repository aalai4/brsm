# Posterior Predictive Draws at New Points

Generates posterior predictive values at arbitrary new predictor points.
Unlike \[predict_surface()\], this function can include residual
uncertainty (sigma), so outputs are draws from the posterior predictive
distribution rather than posterior means alone.

## Usage

``` r
posterior_predict_brsm(
  object,
  factor_names,
  newdata,
  include_residual = TRUE,
  sigma = NULL,
  summary = FALSE,
  probs = c(0.025, 0.5, 0.975),
  draw_subset = NULL,
  max_draws = NULL,
  return_matrix = FALSE,
  output_chunk_size = NULL,
  seed = NULL,
  .sigma_source_draws = NULL
)

# S3 method for class 'brsm_fit'
posterior_predict_brsm(
  object,
  factor_names = NULL,
  newdata,
  include_residual = TRUE,
  sigma = NULL,
  summary = FALSE,
  probs = c(0.025, 0.5, 0.975),
  draw_subset = NULL,
  max_draws = NULL,
  return_matrix = FALSE,
  output_chunk_size = NULL,
  seed = NULL,
  .sigma_source_draws = NULL
)

# S3 method for class 'brmsfit'
posterior_predict_brsm(
  object,
  factor_names,
  newdata,
  include_residual = TRUE,
  sigma = NULL,
  summary = FALSE,
  probs = c(0.025, 0.5, 0.975),
  draw_subset = NULL,
  max_draws = NULL,
  return_matrix = FALSE,
  output_chunk_size = NULL,
  seed = NULL,
  .sigma_source_draws = NULL
)

# Default S3 method
posterior_predict_brsm(
  object,
  factor_names,
  newdata,
  include_residual = TRUE,
  sigma = NULL,
  summary = FALSE,
  probs = c(0.025, 0.5, 0.975),
  draw_subset = NULL,
  max_draws = NULL,
  return_matrix = FALSE,
  output_chunk_size = NULL,
  seed = NULL,
  ...
)
```

## Arguments

- object:

  A `brsm_fit` object, `brmsfit` object, or data frame of posterior
  draws.

- factor_names:

  Character vector of factor names. For `brsm_fit`, defaults to metadata
  when omitted.

- newdata:

  Data frame of predictor values at which to predict.

- include_residual:

  Logical; if `TRUE` (default), include residual noise using posterior
  `sigma` draws.

- sigma:

  Optional residual SD specification. If `NULL`, the function attempts
  to infer `sigma` from posterior draw columns. Can be a scalar or one
  value per posterior draw.

- summary:

  Logical; if `TRUE`, return posterior summaries per row.

- probs:

  Probabilities used for quantile summaries when `summary = TRUE`.

- draw_subset:

  Optional logical or numeric draw subset.

- max_draws:

  Optional cap on number of draws.

- return_matrix:

  Logical; if `TRUE` and `summary = FALSE`, return a draw-by-point
  matrix instead of long-format data.

- output_chunk_size:

  Optional chunk size for long-format assembly when `summary = FALSE`
  and `return_matrix = FALSE`.

- seed:

  Optional random seed for predictive noise draws.

- .sigma_source_draws:

  Internal. Raw posterior draw data frame used to resolve sigma when
  called from `brsm_fit` or `brmsfit` methods. Users should not set this
  directly.

- ...:

  Additional arguments (currently unused).

## Value

A prediction object as matrix, long-format data frame, or summary data
frame depending on options.

## Details

Supports S3 dispatch for `brsm_fit`, `brmsfit`, and data frames of
posterior draws.
