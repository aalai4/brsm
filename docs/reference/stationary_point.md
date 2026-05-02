# Identify Stationary Points of Response Surface

Computes stationary (critical) points of the response surface from
posterior draws, including computation of the eigenvalues of the Hessian
matrix for classification.

## Usage

``` r
stationary_point(
  object,
  factor_names = NULL,
  kappa_thresh = 1e+10,
  diagnostics = c("none", "basic", "full"),
  auto_guidance = TRUE,
  sensitivity_thresholds = c(1e+08, 1e+10, 1e+12)
)

# S3 method for class 'brsm_fit'
stationary_point(
  object,
  factor_names = NULL,
  kappa_thresh = 1e+10,
  diagnostics = c("none", "basic", "full"),
  auto_guidance = TRUE,
  sensitivity_thresholds = c(1e+08, 1e+10, 1e+12)
)

# Default S3 method
stationary_point(
  object,
  factor_names = NULL,
  kappa_thresh = 1e+10,
  diagnostics = c("none", "basic", "full"),
  auto_guidance = TRUE,
  sensitivity_thresholds = c(1e+08, 1e+10, 1e+12)
)
```

## Arguments

- object:

  A \`brsm_fit\` object, \`brmsfit\` object, or data frame of posterior
  draws with Bayesian coefficient columns.

- factor_names:

  Character vector of factor names (if object is data frame).

- kappa_thresh:

  Threshold for computing condition numbers of Hessian.

- diagnostics:

  Diagnostic reporting level. One of `"none"` (default), `"basic"`, or
  `"full"`. When not `"none"`, a diagnostics list is attached as an
  attribute on the returned data frame.

- auto_guidance:

  Logical; if `TRUE` (default), emit guidance warnings when exclusion
  rates are high or threshold sensitivity is unstable.

- sensitivity_thresholds:

  Numeric vector of condition-number thresholds used when
  `diagnostics = "full"`. Defaults to `c(1e8, 1e10, 1e12)`.

## Examples

``` r
if (FALSE) { # \dontrun{
# Fit a model first
fit <- fit_brsm(
  data = my_data,
  response = "y",
  factor_names = c("x1", "x2"),
  chains = 2, iter = 1000, seed = 42
)

# Compute posterior stationary points (x* = -1/2 B^{-1} b)
sp <- stationary_point(fit)
head(sp)

# Posterior mean location of the optimum
colMeans(sp, na.rm = TRUE)

# With basic diagnostics to see exclusion rate
sp <- stationary_point(fit, diagnostics = "basic")
attr(sp, "diagnostics")
} # }
```
