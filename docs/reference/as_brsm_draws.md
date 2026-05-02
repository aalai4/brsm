# Convert Model Output to brsm Draws Format

Converts posterior draws to the brsm-compatible format with standardized
column names (b_Intercept, b_x1, b_I(x1^2), etc.).

## Usage

``` r
as_brsm_draws(object, factor_names, ...)

# S3 method for class 'brsm_fit'
as_brsm_draws(object, factor_names = NULL, ...)

# S3 method for class 'brmsfit'
as_brsm_draws(object, factor_names = NULL, ...)

# S3 method for class 'data.frame'
as_brsm_draws(
  object,
  factor_names,
  require_quadratic = TRUE,
  require_interactions = TRUE,
  ...
)

# Default S3 method
as_brsm_draws(object, factor_names, ...)
```

## Arguments

- object:

  A \`brsm_fit\` object, \`brmsfit\` object, or a data frame of
  posterior draws/coefficients.

- factor_names:

  Character vector of factor names (e.g., c("x1", "x2")).

- ...:

  Additional arguments (currently unused).

- require_quadratic:

  Logical; if `TRUE` (default), require quadratic terms (`b_I(x^2)`
  columns) to be present in the draw data frame.

- require_interactions:

  Logical; if `TRUE` (default), require interaction term columns to be
  present in the draw data frame.

## Value

A data frame with brsm-compatible column names: - b_Intercept:
intercept - b\_\<factor\>: linear term for each factor -
b_I(\<factor\>^2): quadratic term for each factor -
b\_\<factor\>:\<factor\>: interaction terms

## Examples

``` r
# From posterior draws in a data frame
draws_raw <- data.frame(
  b_Intercept = rnorm(100, 5, 0.2),
  b_x1 = rnorm(100, 2, 0.1),
  b_x2 = rnorm(100, 4, 0.1),
  "b_I(x1^2)" = rnorm(100, -1, 0.05),
  "b_I(x2^2)" = rnorm(100, -2, 0.05),
  "b_x1:x2" = rnorm(100, 0.5, 0.05),
  check.names = FALSE
)
draws <- as_brsm_draws(draws_raw, factor_names = c("x1", "x2"))
head(draws)
#>   b_Intercept     b_x1     b_x2  b_I(x1^2) b_I(x2^2)   b_x1:x2
#> 1    5.051063 1.921457 4.136046 -1.0532232 -2.042262 0.5331089
#> 2    4.512547 1.894326 3.992914 -0.9461442 -2.048079 0.5145565
#> 3    4.998886 1.920446 3.972785 -0.9409212 -1.949125 0.5098979
#> 4    5.124311 1.824372 3.755332 -0.9900804 -2.074803 0.4398217
#> 5    5.229682 1.930946 4.006549 -1.0200203 -2.059241 0.4980091
#> 6    4.635636 1.944146 3.890149 -0.9691923 -1.968488 0.5343491
```
