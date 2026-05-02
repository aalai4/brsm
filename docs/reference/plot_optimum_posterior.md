# Plot Posterior Distribution of the Optimum

Plot Posterior Distribution of the Optimum

## Usage

``` r
plot_optimum_posterior(
  draws,
  factor_names,
  stationary_draws = NULL,
  bins = 60,
  levels = 10,
  alpha = 0.9,
  show_points = FALSE,
  point_alpha = 0.2,
  seed = NULL
)
```

## Arguments

- draws:

  Validated draw data frame.

- factor_names:

  Character vector of factor names (exactly 2).

- stationary_draws:

  Optional precomputed stationary draws.

- bins:

  Number of density bins.

- levels:

  Number of density contour levels.

- alpha:

  Transparency for density fill.

- show_points:

  Logical; if TRUE, overlay posterior sample points.

- point_alpha:

  Transparency for point overlay.

- seed:

  Random seed.

## Value

A ggplot2 object.

## Details

Density fills use the viridis "turbo" palette, which is perceptually
uniform and colorblind-safe. Density contour lines are drawn in black
for contrast. The posterior mean optimum is marked with a red cross
symbol. This plot is only available for optimization_geometry mode (when
quadratic terms are present).
