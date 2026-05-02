# brsm: Bayesian Response Surface Methods

------------------------------------------------------------------------

editor_options: markdown: wrap: 72 —

Tools for analyzing posterior distributions of quadratic response
surfaces from Bayesian model fits. The package focuses on Bayesian RSM
workflows: fitting quadratic models, posterior surface prediction,
stationary-point diagnostics, ridge/ascent analysis, and model-checking
utilities.

## Features

- **Bayesian quadratic model fitting** with
  [`fit_brsm()`](https://aalai4.github.io/brsm/reference/fit_brsm.md)
- **Posterior surface prediction** with
  [`predict_surface()`](https://aalai4.github.io/brsm/reference/predict_surface.md)
  and
  [`posterior_predict_brsm()`](https://aalai4.github.io/brsm/reference/posterior_predict_brsm.md)
- **Canonical/stationarity diagnostics** with
  [`canonical_analysis()`](https://aalai4.github.io/brsm/reference/canonical_analysis.md),
  [`stationary_point()`](https://aalai4.github.io/brsm/reference/stationary_point.md),
  and
  [`classify_stationary_point()`](https://aalai4.github.io/brsm/reference/classify_stationary_point.md)
- **Optimization tools** including
  [`posterior_ridge_analysis()`](https://aalai4.github.io/brsm/reference/posterior_ridge_analysis.md),
  [`steepest_ascent()`](https://aalai4.github.io/brsm/reference/steepest_ascent.md),
  [`credible_optimum_region()`](https://aalai4.github.io/brsm/reference/credible_optimum_region.md),
  and
  [`optimize_brsm_multiresponse()`](https://aalai4.github.io/brsm/reference/optimize_brsm_multiresponse.md)
- **Design and priors helpers** via
  [`generate_brsm_design()`](https://aalai4.github.io/brsm/reference/generate_brsm_design.md)
  and
  [`specify_brsm_priors()`](https://aalai4.github.io/brsm/reference/specify_brsm_priors.md)
- **Model adequacy checks** with
  [`loftest_brsm()`](https://aalai4.github.io/brsm/reference/loftest_brsm.md),
  [`check_brsm_fit()`](https://aalai4.github.io/brsm/reference/check_brsm_fit.md),
  and
  [`check_brsm_ppc()`](https://aalai4.github.io/brsm/reference/check_brsm_ppc.md)

## Installation

You can install the development version of brsm from GitHub:

``` r

# install.packages("devtools")
devtools::install_github("aalai4/brsm")
```

## Quick Start

### Basic Workflow

``` r

library(brsm)

# Fit a Bayesian response surface model
fit <- fit_brsm(
  data = my_data,
  response = "y",
  factor_names = c("x1", "x2"),
  chains = 4,
  iter = 2000
)

# Characterize the response surface
stationarity <- fit |> 
  stationary_point() |>
  classify_stationary_point()

ridge <- fit |> posterior_ridge_analysis()

# Predict at new design points including residual uncertainty
new_points <- data.frame(x1 = c(-1, 0, 1), x2 = c(-1, 0, 1))
ppd <- posterior_predict_brsm(fit, newdata = new_points, summary = TRUE)
```

### S3 Method Dispatch

Most analysis functions can take either a `brsm_fit` object or posterior
draws:

``` r

# Recommended: pass brsm_fit directly
sp <- stationary_point(fit)

# Alternate: pass standardized draws explicitly
draws <- as_brsm_draws(fit)
sp2 <- stationary_point(draws, factor_names = c("x1", "x2"))
```

## Core Functions

### Model Fitting

- [`fit_brsm()`](https://aalai4.github.io/brsm/reference/fit_brsm.md):
  Fit Bayesian quadratic response surface model
- [`prepare_brsm_data()`](https://aalai4.github.io/brsm/reference/prepare_brsm_data.md):
  Prepare and center/scale factor variables

### Analysis Functions (S3 Methods)

- [`canonical_analysis()`](https://aalai4.github.io/brsm/reference/canonical_analysis.md):
  Posterior canonical decomposition of the quadratic form
- [`stationary_point()`](https://aalai4.github.io/brsm/reference/stationary_point.md):
  Identify critical points of the response surface
- [`classify_stationary_point()`](https://aalai4.github.io/brsm/reference/classify_stationary_point.md):
  Classify points (maximum, minimum, saddle)
- [`posterior_ridge_analysis()`](https://aalai4.github.io/brsm/reference/posterior_ridge_analysis.md):
  Ridge analysis at specified radii
- [`credible_optimum_region()`](https://aalai4.github.io/brsm/reference/credible_optimum_region.md):
  Compute credible regions around optimum
- [`steepest_ascent()`](https://aalai4.github.io/brsm/reference/steepest_ascent.md):
  Compute steepest ascent path from starting point
- [`posterior_predict_brsm()`](https://aalai4.github.io/brsm/reference/posterior_predict_brsm.md):
  Posterior predictive draws at new points
- [`optimize_brsm_multiresponse()`](https://aalai4.github.io/brsm/reference/optimize_brsm_multiresponse.md):
  Multi-response desirability optimization
- [`loftest_brsm()`](https://aalai4.github.io/brsm/reference/loftest_brsm.md):
  Lack-of-fit test against more complex reference models

### Utilities

- [`as_brsm_draws()`](https://aalai4.github.io/brsm/reference/as_brsm_draws.md):
  Convert various inputs to standardized draws format
- [`compare_brsm_models()`](https://aalai4.github.io/brsm/reference/compare_brsm_models.md):
  Compare multiple models via LOO/WAIC
- [`predict_surface()`](https://aalai4.github.io/brsm/reference/predict_surface.md):
  Generate predictions across response surface grid
- [`surface_grid()`](https://aalai4.github.io/brsm/reference/surface_grid.md):
  Create grid for response surface visualization
- [`decode_brsm_data()`](https://aalai4.github.io/brsm/reference/decode_brsm_data.md):
  Reverse variable transformations to original scale
- [`get_brsm_coding()`](https://aalai4.github.io/brsm/reference/get_brsm_coding.md):
  Extract coding metadata
- [`specify_brsm_priors()`](https://aalai4.github.io/brsm/reference/specify_brsm_priors.md):
  Build `brms` prior specifications
- [`generate_brsm_design()`](https://aalai4.github.io/brsm/reference/generate_brsm_design.md):
  Generate CCD/BBD designs in coded or natural units

## Example: Full Analysis Pipeline

``` r

library(brsm)

# 1. PREPARE DATA
prepared <- prepare_brsm_data(
  data = experimental_data,
  factor_names = c("temperature", "pressure"),
  method = "zscore"  # or "range"
)

# 2. FIT MODELS
fit <- fit_brsm(
  data = prepared,
  response = "yield",
  factor_names = c("temperature", "pressure")
)

# 3. COMPARE MODELS
comparison <- compare_brsm_models(
  models = list(
    quadratic = fit
  ),
  criterion = "loo"
)

# 4. CHARACTERIZE RESPONSE SURFACE
critical <- fit |> stationary_point()
classification <- critical |> classify_stationary_point()
ridge <- fit |> posterior_ridge_analysis()

# 5. GENERATE PREDICTIONS
grid <- surface_grid(
  ranges = list(temperature = c(50, 150), pressure = c(1, 5)),
  n = 20
)
predictions <- predict_surface(
  draw = as_brsm_draws(fit),
  factor_names = c("temperature", "pressure"),
  newdata = grid,
  summary = TRUE
)

# 6. MULTI-RESPONSE OPTIMIZATION
opt <- optimize_brsm_multiresponse(
  object = fit,
  responses = list(
    yield = list(goal = "max")
  ),
  ranges = list(temperature = c(-2, 2), pressure = c(-2, 2))
)

# 7. DECODE TO ORIGINAL SCALE
predictions_original <- decode_brsm_data(
  predictions[, c("temperature", "pressure")],
  coding = get_brsm_coding(fit)
)

# 8. LACK-OF-FIT TEST
lof_test <- loftest_brsm(
  object = fit,
  reference_type = "cubic",
  include_ppc = TRUE
)
```

## Design Philosophy

The package architecture prioritizes:

- **Usability**: Intuitive, pipe-compatible API

- **Bayesian Rigor**: Full posterior uncertainty for all estimates

- **Backward Compatibility**: Existing code continues to work

- **Flexibility**: Extensible through S3 methods - **Documentation**:
  Comprehensive examples and vignettes

## References

- Guo, X., Luh, D. B., & Box, G. E. (2009). Bayesian non-parametric
  modelling for case studies in operations research. *Journal of the
  Royal Statistical Society: Series C*, 58(1), 99-118.

## Documentation

See package help pages and examples in the `man/` directory for current
function-level documentation.

## License

MIT License - see LICENSE file for details

## Contributing

Contributions are welcome! Please submit issues and pull requests to the
GitHub repository.
