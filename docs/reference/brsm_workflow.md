# Run a Complete BRSM Analysis Workflow

Runs a response-surface workflow from model/draw conversion through
prediction, diagnostics, and optional plots.

## Usage

``` r
brsm_workflow(
  object,
  factor_names,
  ranges = NULL,
  probs = c(0.025, 0.5, 0.975),
  grid_n = 40,
  radii = seq(0, 1, length.out = 5),
  start = NULL,
  step_size = 0.1,
  n_steps = 10,
  include_plots = TRUE,
  contour_type = c("mean", "uncertainty", "quantile"),
  plot_vary_factors = NULL,
  plot_fixed = NULL,
  plot_conditioning = c("center", "optimum_mean", "user"),
  plot_slice = NULL,
  plot_pairwise = FALSE,
  seed = NULL,
  fit_mode = c("none", "fit_brsm"),
  response = NULL,
  fit_args = list(),
  steps = c("predictions", "stationary", "classification", "credible_region",
    "steepest_ascent", "ridge", "plots")
)
```

## Arguments

- object:

  Either (1) a `brsm_fit` object from \[fit_brsm()\], (2) a `brmsfit`
  object from `brms::brm()`, (3) a posterior-draw data frame with
  Bayesian coefficient columns, or (4) raw analysis data when
  `fit_mode != "none"`.

- factor_names:

  Character vector of factor names.

- ranges:

  Named list with factor ranges.

- probs:

  Posterior probabilities for prediction summaries.

- grid_n:

  Number of points per factor dimension in the prediction grid.

- radii:

  Numeric vector of radii for ridge analysis.

- start:

  Optional named numeric vector of starting values for steepest ascent.
  Defaults to zeros for each factor.

- step_size:

  Step size for steepest ascent.

- n_steps:

  Number of steps for steepest ascent.

- include_plots:

  Logical; if \`TRUE\`, generate contour/optimum/ridge plots when
  exactly two factors are provided. For non-quadratic draw objects,
  contour plots are returned with a steepest-ascent direction overlay
  and optimization-geometry plots are omitted.

- contour_type:

  Surface type for contour plotting. Passed to
  \[plot_posterior_contours()\].

- plot_vary_factors:

  Optional character vector of length 2 specifying which factors vary in
  contour plots when \`length(factor_names) \> 2\` and \`plot_pairwise =
  FALSE\`. Defaults to the first two factors.

- plot_fixed:

  Optional named numeric vector of fixed values used when
  \`plot_conditioning = "user"\`.

- plot_conditioning:

  Conditioning strategy for non-varied factors in contour plots. One of
  \`"center"\`, \`"optimum_mean"\`, or \`"user"\`.

- plot_slice:

  Optional named list with one element specifying sliced values for a
  conditioning factor, e.g. \`list(x3 = c(-1, 0, 1))\`.

- plot_pairwise:

  Logical; if \`TRUE\`, generate pairwise contour panels for all factor
  pairs.

- seed:

  Optional random seed used for stochastic plotting operations.

- fit_mode:

  One of \`"none"\` or \`"fit_brsm"\`. Use \`"none"\` (default) to
  analyze an existing Bayesian fit/draw object. Use \`"fit_brsm"\` to
  fit from raw data in `object` and then run the workflow.

- response:

  Required when \`fit_mode != "none"\`. Name of the response column in
  raw input data.

- fit_args:

  Named list of additional arguments forwarded to \[fit_brsm()\] when
  \`fit_mode != "none"\`. Useful for custom priors and sampling
  controls.

- steps:

  Character vector selecting which workflow stages to run. Any of
  \`"predictions"\`, \`"stationary"\`, \`"classification"\`,
  \`"credible_region"\`, \`"steepest_ascent"\`, \`"ridge"\`,
  \`"plots"\`. Default runs all stages.

## Value

A named list with components: always \`draws\`, and conditionally
\`fit\`, \`ranges\`, \`grid\`, \`predictions\`, \`surface\`,
\`stationary_points\`, \`classification\`, \`credible_region\`,
\`steepest_ascent\`, \`ridge_analysis\`, and \`plots\` depending on
\`steps\`. When plots are produced, \`plots\` also includes
\`plot_mode\` (\`"optimization_geometry"\` or \`"direction_only"\`) and
\`geometry_available\` (logical).

## Details

By default (`fit_mode = "none"`), input must be Bayesian: `brsm_fit`,
`brmsfit`, or posterior-draw data frame with required coefficient
columns (`b_Intercept`, `b_x1`, `b_I(x1^2)`, interactions, etc.).

If `fit_mode` is `"fit_brsm"`, `object` is treated as raw input data and
fitting is performed before selected workflow steps run.

## Examples

``` r
set.seed(1)
dat <- data.frame(
  x1 = runif(50, -2, 2),
  x2 = runif(50, -2, 2)
)
dat$y <- 5 + 2 * dat$x1 + 4 * dat$x2 - dat$x1^2 - 2 * dat$x2^2 +
  0.5 * dat$x1 * dat$x2 + rnorm(50, sd = 0.5)
draws_raw <- data.frame(
  b_Intercept = rnorm(200, 5, 0.2),
  b_x1 = rnorm(200, 2, 0.1),
  b_x2 = rnorm(200, 4, 0.1),
  "b_I(x1^2)" = rnorm(200, -1, 0.05),
  "b_I(x2^2)" = rnorm(200, -2, 0.05),
  "b_x1:x2" = rnorm(200, 0.5, 0.05),
  check.names = FALSE
)

out <- brsm_workflow(
  object = draws_raw,
  factor_names = c("x1", "x2"),
  ranges = list(x1 = c(-2, 2), x2 = c(-2, 2))
)

names(out)
#>  [1] "draws"             "ranges"            "grid"             
#>  [4] "predictions"       "surface"           "stationary_points"
#>  [7] "classification"    "credible_region"   "steepest_ascent"  
#> [10] "ridge_analysis"    "plots"            
```
