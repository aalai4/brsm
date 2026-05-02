# Compute Credible Region Around Response Surface Optimum

Identifies the credible region around the stationary point of the
response surface at specified probability levels. Useful for uncertainty
quantification around the estimated optimum.

## Usage

``` r
credible_optimum_region(
  object,
  factor_names = NULL,
  probs = c(0.025, 0.975),
  summary = TRUE
)

# S3 method for class 'brsm_fit'
credible_optimum_region(
  object,
  factor_names = NULL,
  probs = c(0.025, 0.975),
  summary = TRUE
)

# Default S3 method
credible_optimum_region(
  object,
  factor_names = NULL,
  probs = c(0.025, 0.975),
  summary = TRUE
)
```

## Arguments

- object:

  A `brsm_fit` object, `brmsfit` object, or data frame of posterior
  draws with Bayesian coefficient columns.

- factor_names:

  Character vector of factor names (if object is data frame).

- probs:

  Probability levels for credible intervals (default: 95% interval).

- summary:

  If `TRUE`, return region boundaries as a summary data frame; if
  `FALSE`, return all posterior draws.

## Value

Data frame with credible region boundaries or posterior draws.
