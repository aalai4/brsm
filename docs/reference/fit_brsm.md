# Fit a Bayesian Quadratic Response Surface Model

Fits a second-order (quadratic) response surface model using
`brms::brm()` and returns a `brsm_fit` object that stores the fitted
model and metadata for downstream brsm analysis.

## Usage

``` r
fit_brsm(
  data,
  response,
  factor_names,
  ranges = NULL,
  prior = NULL,
  prior_profile = c("legacy", "regularized", "adaptive"),
  family = stats::gaussian(),
  chains = 4,
  iter = 2000,
  warmup = floor(iter/2),
  seed = NULL,
  sampling_preset = c("fast", "balanced", "robust"),
  backend = NULL,
  control = NULL,
  model_terms = c("second_order", "first_order", "first_order_twi", "pure_quadratic"),
  coding_policy = c("warn", "error", "ignore"),
  ...
)
```

## Arguments

- data:

  Data frame containing response and factor columns.

- response:

  Name of the response variable.

- factor_names:

  Character vector of factor names.

- ranges:

  Optional named list of factor ranges. If `NULL`, ranges are inferred
  from `data`.

- prior:

  Optional prior specification passed to `brms::brm()`.

- prior_profile:

  Prior profile to use when `prior = NULL`. One of `"legacy"` (flat,
  wide normal priors), `"regularized"` (moderate shrinkage), or
  `"adaptive"` (data-scaled priors). Default is `"legacy"`.

- family:

  Model family passed to `brms::brm()`.

- chains:

  Number of MCMC chains.

- iter:

  Number of iterations per chain.

- warmup:

  Number of warmup iterations per chain.

- seed:

  Optional random seed.

- sampling_preset:

  Sampling profile for NUTS tuning. One of `"fast"`, `"balanced"`, or
  `"robust"`.

- backend:

  Optional backend passed to `brms::brm()` (e.g., `"rstan"`,
  `"cmdstanr"`). If `NULL`, uses brms defaults.

- control:

  Optional named list of NUTS control arguments. Values here override
  defaults from `sampling_preset`.

- model_terms:

  Polynomial term specification. One of `"second_order"` (default;
  linear + two-way interactions + pure quadratic terms), `"first_order"`
  (linear only), `"first_order_twi"` (linear + two-way interactions), or
  `"pure_quadratic"` (linear + pure quadratic terms).

- coding_policy:

  How to handle missing coding metadata from \[prepare_brsm_data()\].
  One of `"warn"` (default), `"error"`, or `"ignore"`.

- ...:

  Additional arguments passed to `brms::brm()`.

## Value

An object of class `brsm_fit` with elements: `fit`, `formula`,
`response`, `factor_names`, `ranges`, `model_terms`, and `call`.

## Examples

``` r
if (FALSE) { # \dontrun{
dat <- data.frame(
  x1 = runif(100, -2, 2),
  x2 = runif(100, -2, 2)
)
dat$y <- 5 + 2 * dat$x1 + 4 * dat$x2 - dat$x1^2 - 2 * dat$x2^2 +
  0.5 * dat$x1 * dat$x2 + rnorm(100, sd = 0.5)

fit <- fit_brsm(
  data = dat,
  response = "y",
  factor_names = c("x1", "x2"),
  sampling_preset = "balanced",
  control = list(adapt_delta = 0.95),
  chains = 2,
  iter = 1000,
  warmup = 500,
  seed = 123
)

out <- brsm_workflow(
  object = fit,
  factor_names = fit$factor_names,
  ranges = fit$ranges
)
} # }
```
