# Compute Steepest Ascent Path on Response Surface

Computes the path of steepest ascent (or descent) on the response
surface starting from a specified point. Useful for sequential design of
experiments to locate the vicinity of the optimum.

## Usage

``` r
steepest_ascent(
  object,
  factor_names = NULL,
  start = NULL,
  step_size = 1,
  n_steps = 10,
  tol = 1e-08,
  return_mean_path = TRUE,
  drop_null_paths = TRUE
)

# S3 method for class 'brsm_fit'
steepest_ascent(
  object,
  factor_names = NULL,
  start = NULL,
  step_size = 1,
  n_steps = 10,
  tol = 1e-08,
  return_mean_path = TRUE,
  drop_null_paths = TRUE
)

# Default S3 method
steepest_ascent(
  object,
  factor_names = NULL,
  start = NULL,
  step_size = 1,
  n_steps = 10,
  tol = 1e-08,
  return_mean_path = TRUE,
  drop_null_paths = TRUE
)
```

## Arguments

- object:

  A \`brsm_fit\` object, \`brmsfit\` object, or data frame of posterior
  draws with Bayesian coefficient columns.

- factor_names:

  Character vector of factor names (if object is data frame).

- start:

  Starting point for ascent path (default: center of coded region).

- step_size:

  Step size along gradient direction.

- n_steps:

  Number of steps to compute.

- tol:

  Numerical tolerance for gradient computations.

- return_mean_path:

  If TRUE, return mean path; if FALSE, return all draws.

- drop_null_paths:

  Drop paths that don't increase response.

## Value

Data frame with steepest ascent path coordinates and response values.
