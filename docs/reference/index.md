# Package index

## Package Overview

- [`brsm-package`](https://aalai4.github.io/brsm/reference/brsm-package.md)
  [`brsm`](https://aalai4.github.io/brsm/reference/brsm-package.md) :
  brsm: Bayesian Response Surface Methods

## Model Fitting

Fit Bayesian second-order response surface models.

- [`fit_brsm()`](https://aalai4.github.io/brsm/reference/fit_brsm.md) :
  Fit a Bayesian Quadratic Response Surface Model
- [`specify_brsm_priors()`](https://aalai4.github.io/brsm/reference/specify_brsm_priors.md)
  : Build Prior Specifications for BRSM Models
- [`check_brsm_fit()`](https://aalai4.github.io/brsm/reference/check_brsm_fit.md)
  : Check MCMC Diagnostics for a brsm Fit
- [`check_brsm_ppc()`](https://aalai4.github.io/brsm/reference/check_brsm_ppc.md)
  : Posterior Predictive Check Summary for a brsm Fit

## Posterior Analysis

Extract and analyse the fitted posterior distribution.

- [`as_brsm_draws()`](https://aalai4.github.io/brsm/reference/as_brsm_draws.md)
  : Convert Model Output to brsm Draws Format
- [`stationary_point()`](https://aalai4.github.io/brsm/reference/stationary_point.md)
  : Identify Stationary Points of Response Surface
- [`hessian_quadratic()`](https://aalai4.github.io/brsm/reference/hessian_quadratic.md)
  : Compute the Posterior Hessian of a Quadratic Response Surface
- [`gradient_quadratic()`](https://aalai4.github.io/brsm/reference/gradient_quadratic.md)
  : Compute the Posterior Gradient of a Quadratic Response Surface
- [`canonical_analysis()`](https://aalai4.github.io/brsm/reference/canonical_analysis.md)
  : Canonical Analysis of Bayesian Response Surface
- [`classify_stationary_point()`](https://aalai4.github.io/brsm/reference/classify_stationary_point.md)
  : Classify Stationary Points as Maximum, Minimum, or Saddle
- [`summarize_brsm_stability()`](https://aalai4.github.io/brsm/reference/summarize_brsm_stability.md)
  : Summarize BRSM Stability Diagnostics
- [`decode_stationary_status()`](https://aalai4.github.io/brsm/reference/decode_stationary_status.md)
  : Decode Stationary Solver Status Codes

## Uncertainty Quantification

Credible regions, ridge analysis, and steepest ascent.

- [`credible_optimum_region()`](https://aalai4.github.io/brsm/reference/credible_optimum_region.md)
  : Compute Credible Region Around Response Surface Optimum
- [`posterior_ridge_analysis()`](https://aalai4.github.io/brsm/reference/posterior_ridge_analysis.md)
  : Analyze Response Surface Along Ridges
- [`steepest_ascent()`](https://aalai4.github.io/brsm/reference/steepest_ascent.md)
  : Compute Steepest Ascent Path on Response Surface
- [`posterior_predict_brsm()`](https://aalai4.github.io/brsm/reference/posterior_predict_brsm.md)
  : Posterior Predictive Draws at New Points

## Surface Prediction

Predict the response surface over factor grids.

- [`predict_surface()`](https://aalai4.github.io/brsm/reference/predict_surface.md)
  : Predict Response Surface from Posterior Draws or brsm_fit Objects
- [`surface_grid()`](https://aalai4.github.io/brsm/reference/surface_grid.md)
  : Generate a Regular Grid over Factor Ranges
- [`plot_posterior_contours()`](https://aalai4.github.io/brsm/reference/plot_posterior_contours.md)
  : Plot Posterior Contours

## Design and Data Helpers

Generate experimental designs and manage coded factor levels.

- [`generate_brsm_design()`](https://aalai4.github.io/brsm/reference/generate_brsm_design.md)
  : Generate Standard Response Surface Designs
- [`prepare_brsm_data()`](https://aalai4.github.io/brsm/reference/prepare_brsm_data.md)
  : Prepare Factor Columns for BRSM Modeling
- [`decode_brsm_data()`](https://aalai4.github.io/brsm/reference/decode_brsm_data.md)
  : Decode Coded Factor Columns Back to Original Scale
- [`get_brsm_coding()`](https://aalai4.github.io/brsm/reference/get_brsm_coding.md)
  : Extract Stored brsm Coding Metadata

## Workflow and Utilities

High-level workflow runner and model comparison tools.

- [`brsm_workflow()`](https://aalai4.github.io/brsm/reference/brsm_workflow.md)
  : Run a Complete BRSM Analysis Workflow
- [`compare_brsm_models()`](https://aalai4.github.io/brsm/reference/compare_brsm_models.md)
  : Compare Bayesian Models for brsm
- [`loftest_brsm()`](https://aalai4.github.io/brsm/reference/loftest_brsm.md)
  : Bayesian Lack-of-Fit Convenience Wrapper for brsm
- [`optimize_brsm_multiresponse()`](https://aalai4.github.io/brsm/reference/optimize_brsm_multiresponse.md)
  : Multi-Response Optimization via Posterior Desirability
