---

editor_options: 
  markdown: 
    wrap: 72
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

# brsm: Bayesian Response Surface Methods

<!-- badges: start -->

<!-- badges: end -->

Tools for analyzing posterior distributions of quadratic response surfaces from Bayesian model fits. The package focuses on Bayesian RSM workflows: fitting quadratic models, posterior surface prediction, stationary-point diagnostics, ridge/ascent analysis, and model-checking utilities.

## Features

- **Bayesian quadratic model fitting** with `fit_brsm()`
- **Posterior surface prediction** with `predict_surface()` and `posterior_predict_brsm()`
- **Canonical/stationarity diagnostics** with `canonical_analysis()`, `stationary_point()`, and `classify_stationary_point()`
- **Optimization tools** including `posterior_ridge_analysis()`, `steepest_ascent()`, `credible_optimum_region()`, and `optimize_brsm_multiresponse()`
- **Design and priors helpers** via `generate_brsm_design()` and `specify_brsm_priors()`
- **Model adequacy checks** with `loftest_brsm()`, `check_brsm_fit()`, and `check_brsm_ppc()`

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

Most analysis functions can take either a `brsm_fit` object or posterior draws:

``` r
# Recommended: pass brsm_fit directly
sp <- stationary_point(fit)

# Alternate: pass standardized draws explicitly
draws <- as_brsm_draws(fit)
sp2 <- stationary_point(draws, factor_names = c("x1", "x2"))
```

## Core Functions

### Model Fitting

- `fit_brsm()`: Fit Bayesian quadratic response surface model
- `prepare_brsm_data()`: Prepare and center/scale factor variables

### Analysis Functions (S3 Methods)

- `canonical_analysis()`: Posterior canonical decomposition of the quadratic form
- `stationary_point()`: Identify critical points of the response surface
- `classify_stationary_point()`: Classify points (maximum, minimum, saddle)
- `posterior_ridge_analysis()`: Ridge analysis at specified radii
- `credible_optimum_region()`: Compute credible regions around optimum
- `steepest_ascent()`: Compute steepest ascent path from starting point
- `posterior_predict_brsm()`: Posterior predictive draws at new points
- `optimize_brsm_multiresponse()`: Multi-response desirability optimization
- `loftest_brsm()`: Lack-of-fit test against more complex reference models

### Utilities

- `as_brsm_draws()`: Convert various inputs to standardized draws format
- `compare_brsm_models()`: Compare multiple models via LOO/WAIC
- `predict_surface()`: Generate predictions across response surface grid
- `surface_grid()`: Create grid for response surface visualization
- `decode_brsm_data()`: Reverse variable transformations to original scale
- `get_brsm_coding()`: Extract coding metadata
- `specify_brsm_priors()`: Build `brms` prior specifications
- `generate_brsm_design()`: Generate CCD/BBD designs in coded or natural units

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

- **Flexibility**: Extensible through S3 methods - **Documentation**: Comprehensive examples and vignettes

## References

- Guo, X., Luh, D. B., & Box, G. E. (2009). Bayesian non-parametric modelling for case studies in operations research. *Journal of the Royal Statistical Society: Series C*, 58(1), 99-118.

## Documentation

See package help pages and examples in the `man/` directory for current function-level documentation.

## License

MIT License - see LICENSE file for details

## Contributing

Contributions are welcome! Please submit issues and pull requests to the GitHub repository.


## For the class presentation:
https://uofnebraska-my.sharepoint.com/:p:/g/personal/65886697_nebraska_edu/IQCJ4ubigj1KRI1wfz63fcn_ATTItva7u-jermFwAu3u454?e=lt0fUH

# Data Prep ----------------------------------------------------------
library(brsm)
library(brms)   # required for fitting

x=scan(
  text="1 1 1.000 1.000 1.000 6.2 30.5 20.0 0.547 5.06 36.35
 2 1 1.000 -1.000 -1.000 6.2 25.5 5.0 0.643 6.05 21.30
 3 1 -1.000 1.000 -1.000 3.8 30.5 5.0 0.112 1.19 15.96
 4 1 -1.000 -1.000 1.000 3.8 25.5 20.0 0.272 3.96 43.60
 5 1 0.000 0.000 0.000 5.0 28.0 12.5 0.442 3.52 32.74
 6 1 0.000 0.000 0.000 5.0 28.0 12.5 0.108 2.02 19.91
 7 2 1.000 1.000 -1.000 6.2 30.5 5.0 0.752 7.50 79.16
 8 2 1.000 -1.000 1.000 6.2 25.5 20.0 0.615 5.12 23.06
 9 2 -1.000 1.000 1.000 3.8 30.5 20.0 0.063 0.77 18.57
 10 2 -1.000 -1.000 -1.000 3.8 25.5 5.0 0.144 1.25 11.08
 11 2 0.000 0.000 0.000 5.0 28.0 12.5 0.341 3.96 25.29
 12 2 0.000 0.000 0.000 5.0 28.0 12.5 0.315 3.20 46.93
 13 0 -1.633 0.000 0.000 3.0 28.0 12.5 0.637 5.89 59.65
 14 0 1.633 0.000 0.000 7.0 28.0 12.5 1.029 9.98 89.93
 15 0 0.000 -1.633 0.000 5.0 24.0 12.5 0.248 1.85 25.58
 16 0 0.000 1.633 0.000 5.0 32.0 12.5 0.008 0.37 4.21
 17 0 0.000 0.000 -1.633 5.0 28.0 0.0 0.024 0.70 3.72
 18 0 0.000 0.000 1.633 5.0 28.0 25.0 0.037 0.38 8.22
 19 0 0.000 0.000 0.000 5.0 28.0 12.5 0.638 5.74 65.31
 20 0 0.000 0.000 0.000 5.0 28.0 12.5 0.375 3.90 42.27 ",
  what=list(o=0, bl=0, xp=0, xt=0, xg=0, p=0, t=0, gly=0, y=0, spy=0, spa=0),nlines=21)

#obs block xpH xtemp xgly pH temp glycerol yield sp.yld sp.act ;
d=data.frame(x)
attach(d)
names(d)
d$block=factor(bl)
dat <- d[, c("xp", "xt", "xg", "y")]

dat_coded <- prepare_brsm_data(
  data = dat,
  factor_names = c("xp", "xt", "xg"),
  method = "identity"
)

fit <- fit_brsm(
  data = dat_coded,
  response = "y",
  factor_names = c("xp", "xt", "xg"),
  chains = 2,
  iter = 2000,
  warmup = 1000,
  seed = 123,
  sampling_preset = "balanced",
  coding_policy = "ignore", # because already coded
  refresh = 0,
  silent = 2
)

print(fit)
print(summary(fit))
print(check_brsm_fit(fit))
print(check_brsm_ppc(fit))

# Stationary Point Analysis -----------------------------------------------
stat_draws <- stationary_point(object = fit, factor_names = c("xp", "xt", "xg"))

# Summary: posterior mean and credible interval
cat("Posterior mean stationary point:\n")
print(colMeans(stat_draws))

cat("\n95% credible interval (pointwise):\n")
print(apply(stat_draws, 2, quantile, probs = c(0.025, 0.975)))

# Visualization and EDA ---------------------------------------------------
draws_df <- as_brsm_draws(fit, factor_names = c("xp", "xt", "xg")) # Build canonical posterior draws once

p_pairs <- brsm:::plot_posterior_contours(
  draws = draws_df,
  factor_names = c("xp", "xt", "xg"),
  ranges = fit$ranges,
  bins = 12,
  pairwise = TRUE
)
print(p_pairs)
