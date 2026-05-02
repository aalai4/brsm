# Build Prior Specifications for BRSM Models

Builds `brms` priors for quadratic response-surface models used by
\[fit_brsm()\]. Priors are generated for intercept, linear, interaction,
quadratic, and residual sigma terms according to `model_terms`.

## Usage

``` r
specify_brsm_priors(
  factor_names,
  model_terms = c("second_order", "first_order", "first_order_twi", "pure_quadratic"),
  prior_profile = c("legacy", "regularized", "adaptive"),
  coefficient_family = c("normal", "student_t"),
  intercept_sd = 5,
  linear_sd = 2,
  interaction_sd = 1,
  quadratic_sd = 1,
  sigma_scale = 2.5,
  student_df = 3,
  include_intercept = TRUE,
  include_sigma = TRUE,
  autoscale = FALSE,
  data = NULL,
  response = NULL
)
```

## Arguments

- factor_names:

  Character vector of factor names.

- model_terms:

  Polynomial term specification. One of `"second_order"`,
  `"first_order"`, `"first_order_twi"`, or `"pure_quadratic"`.

- prior_profile:

  Prior profile controlling default scale parameters. One of `"legacy"`
  (wide flat priors), `"regularized"` (moderately shrinking), or
  `"adaptive"` (data-scaled). User-supplied scale arguments always
  override profile defaults.

- coefficient_family:

  Prior family for intercept and slope terms. One of `"normal"` or
  `"student_t"`.

- intercept_sd:

  Base scale for intercept prior.

- linear_sd:

  Base scale for linear term priors.

- interaction_sd:

  Base scale for interaction term priors.

- quadratic_sd:

  Base scale for quadratic term priors.

- sigma_scale:

  Base scale for residual `sigma` prior.

- student_df:

  Degrees of freedom when `coefficient_family = "student_t"`.

- include_intercept:

  Logical; include intercept prior.

- include_sigma:

  Logical; include residual sigma prior.

- autoscale:

  Logical; if `TRUE` and `data`/`response` are provided, scales all
  prior standard deviations by `sd(data[[response]])`.

- data:

  Optional data frame used for autoscaling.

- response:

  Optional response column name used for autoscaling.

## Value

A `brmsprior` object that can be passed directly to \[fit_brsm()\] or
`brms::brm()` as the `prior` argument.
