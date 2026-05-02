# Predict Response Surface from Posterior Draws or brsm_fit Objects

Supports S3 dispatch for \`brsm_fit\`, \`brmsfit\`, and draw data
frames. For \`brsm_fit\` objects, \`factor_names\` are inferred from
metadata when omitted.

## Usage

``` r
predict_surface(
  draws,
  factor_names,
  newdata,
  summary = FALSE,
  probs = c(0.025, 0.5, 0.975),
  draw_subset = NULL,
  max_draws = NULL,
  return_matrix = FALSE,
  output_chunk_size = NULL
)

# S3 method for class 'brsm_fit'
predict_surface(
  draws,
  factor_names = NULL,
  newdata,
  summary = FALSE,
  probs = c(0.025, 0.5, 0.975),
  draw_subset = NULL,
  max_draws = NULL,
  return_matrix = FALSE,
  output_chunk_size = NULL
)

# S3 method for class 'brmsfit'
predict_surface(
  draws,
  factor_names,
  newdata,
  summary = FALSE,
  probs = c(0.025, 0.5, 0.975),
  draw_subset = NULL,
  max_draws = NULL,
  return_matrix = FALSE,
  output_chunk_size = NULL
)

# Default S3 method
predict_surface(
  draws,
  factor_names,
  newdata,
  summary = FALSE,
  probs = c(0.025, 0.5, 0.975),
  draw_subset = NULL,
  max_draws = NULL,
  return_matrix = FALSE,
  output_chunk_size = NULL
)
```

## Arguments

- draws:

  A \`brsm_fit\` object, \`brmsfit\` object, or data frame of draws.

- factor_names:

  Character vector of factor names.

- newdata:

  Data frame of predictor values.

- summary:

  Logical; if \`TRUE\`, return posterior summaries per row.

- probs:

  Probabilities used for quantile summaries when \`summary=TRUE\`.

- draw_subset:

  Optional logical or numeric draw subset.

- max_draws:

  Optional cap on number of draws.

- return_matrix:

  Logical; if \`TRUE\` and \`summary=FALSE\`, return draw-by-point
  matrix instead of long-format data frame.

- output_chunk_size:

  Optional chunk size for long-format assembly.

## Value

A prediction object as matrix, long-format data frame, or summary data
frame depending on options.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- fit_brsm(data = my_data, response = "y",
                factor_names = c("x1", "x2"), chains = 2, seed = 1)

# Define a grid of predictor values
grid <- expand.grid(x1 = seq(-1, 1, length.out = 10),
                    x2 = seq(-1, 1, length.out = 10))

# Posterior mean and 95% credible interval at each grid point
preds <- predict_surface(fit, newdata = grid, summary = TRUE)
head(preds)

# Draw-level predictions (long format)
preds_long <- predict_surface(fit, newdata = grid, summary = FALSE)
} # }
```
