# Plot Posterior Contours

Plot Posterior Contours

## Usage

``` r
plot_posterior_contours(
  draws,
  factor_names,
  ranges,
  n = 60,
  type = c("mean", "uncertainty", "quantile"),
  probs = c(0.025, 0.975),
  quantile = 0.5,
  bins = 10,
  vary_factors = NULL,
  fixed = NULL,
  conditioning = c("center", "optimum_mean", "user"),
  slice = NULL,
  pairwise = FALSE,
  overlay_stationary = FALSE,
  stationary_draws = NULL,
  overlay_type = c("mean", "posterior"),
  overlay_alpha = 0.1,
  overlay_max_draws = 2000,
  seed = NULL
)
```

## Arguments

- draws:

  Validated draw data frame.

- factor_names:

  Character vector of factor names.

- ranges:

  Named list of factor ranges.

- n:

  Grid resolution.

- type:

  Surface type ("mean", "uncertainty", "quantile").

- probs:

  Probability quantiles for uncertainty.

- quantile:

  Quantile value if type="quantile".

- bins:

  Number of contour bins.

- vary_factors:

  Optional 2-factor pair for conditional plotting.

- fixed:

  Optional named numeric vector for fixed factor values.

- conditioning:

  Strategy for fixed factors.

- slice:

  Optional named list for sliced factor values.

- pairwise:

  Logical; if TRUE, generate all factor-pair panels.

- overlay_stationary:

  Logical; if TRUE, overlay stationary points.

- stationary_draws:

  Optional precomputed stationary draws.

- overlay_type:

  Type of stationary overlay ("mean" or "posterior").

- overlay_alpha:

  Transparency for posterior draw overlay.

- overlay_max_draws:

  Cap on posterior draws for overlay.

- seed:

  Random seed.

## Value

A ggplot2 object.

## Details

Contour fills use the viridis "turbo" palette, which is perceptually
uniform and colorblind-safe. Stationary points are overlaid in red
(single-pair, non-sliced contours only). For 3+ factors, use
\`pairwise=TRUE\` for all pairs or \`slice\` for conditional slices.
