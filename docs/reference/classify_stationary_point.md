# Classify Stationary Points as Maximum, Minimum, or Saddle

Classifies stationary points based on the eigenvalues of the Hessian
matrix. Uses posterior distribution to compute probabilities of each
classification.

## Usage

``` r
classify_stationary_point(
  object,
  factor_names = NULL,
  tol = 1e-08,
  return_kappa = FALSE
)

# S3 method for class 'brsm_fit'
classify_stationary_point(
  object,
  factor_names = NULL,
  tol = 1e-08,
  return_kappa = FALSE
)

# Default S3 method
classify_stationary_point(
  object,
  factor_names = NULL,
  tol = 1e-08,
  return_kappa = FALSE
)
```

## Arguments

- object:

  A \`brsm_fit\` object, \`brmsfit\` object, or data frame of posterior
  draws with Bayesian coefficient columns.

- factor_names:

  Character vector of factor names (if object is data frame).

- tol:

  Numerical tolerance for eigenvalue comparisons.

- return_kappa:

  Logical; if `TRUE`, attach a `kappa` column to the output containing
  the LU-diagonal condition-number proxy of the Hessian for each draw.
  Useful for identifying near-singular draws without altering
  classification behaviour. Default `FALSE`.

## Value

A data frame with columns `draw` and `classification` (a factor with
levels `"maximum"`, `"minimum"`, `"saddle"`, `"indeterminate"`). When
`return_kappa = TRUE` an additional numeric column `kappa` is included.
