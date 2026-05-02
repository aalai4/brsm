# Analyze Response Surface Along Ridges

Performs ridge analysis by locating maximum and minimum points along
specified radii from the center of the experimental region. Computes
posterior distributions of ridge values.

## Usage

``` r
posterior_ridge_analysis(
  object,
  factor_names = NULL,
  radii = seq(0, 3, length.out = 10),
  tol = 1e-06,
  summary = TRUE
)

# S3 method for class 'brsm_fit'
posterior_ridge_analysis(
  object,
  factor_names = NULL,
  radii = seq(0, 3, length.out = 10),
  tol = 1e-06,
  summary = TRUE
)

# Default S3 method
posterior_ridge_analysis(
  object,
  factor_names = NULL,
  radii = seq(0, 3, length.out = 10),
  tol = 1e-06,
  summary = TRUE
)
```

## Arguments

- object:

  A \`brsm_fit\` object, \`brmsfit\` object, or data frame of posterior
  draws with Bayesian coefficient columns.

- factor_names:

  Character vector of factor names (if object is data frame).

- radii:

  Sequence of radii at which to compute ridge values.

- tol:

  Numerical tolerance for ridge computations.

- summary:

  If TRUE, return summarized ridge values; if FALSE, return draws.

## Value

Data frame with posterior ridge values or summary statistics.
