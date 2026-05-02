# Bayesian Lack-of-Fit Convenience Wrapper for brsm

Runs a Bayesian lack-of-fit style check by comparing a baseline model
against a richer reference model using LOO/WAIC, with optional
side-by-side posterior predictive summaries.

## Usage

``` r
loftest_brsm(
  object,
  reference_model = NULL,
  reference_type = c("cubic", "extended"),
  criterion = c("loo", "waic"),
  loo_moment_match = FALSE,
  loo_reloo = FALSE,
  loo_k_threshold = 0.7,
  loo_auto_moment_match = TRUE,
  include_ppc = FALSE,
  ppc_ndraws = 200,
  ppc_probs = c(0.025, 0.975),
  seed = NULL,
  data = NULL,
  response = NULL,
  factor_names = NULL,
  prior = NULL,
  family = NULL,
  chains = 4,
  iter = 2000,
  warmup = floor(iter/2),
  sampling_preset = c("fast", "balanced", "robust"),
  backend = NULL,
  control = NULL,
  ...
)
```

## Arguments

- object:

  Baseline model as a `brsm_fit` or `brmsfit` object.

- reference_model:

  Optional richer reference model as `brsm_fit` or `brmsfit`. If `NULL`,
  the function attempts to fit one.

- reference_type:

  Type of richer reference formula used when `reference_model = NULL`.
  One of `"cubic"` (adds cubic terms) or `"extended"` (adds cubic and
  quadratic-by-linear terms).

- criterion:

  Information criterion for model comparison: `"loo"` (default) or
  `"waic"`.

- loo_moment_match:

  Logical; when `criterion = "loo"`, pass `moment_match` to
  \[brms::loo()\] via \[compare_brsm_models()\].

- loo_reloo:

  Logical; when `criterion = "loo"`, pass `reloo` to \[brms::loo()\] via
  \[compare_brsm_models()\].

- loo_k_threshold:

  Pareto-k threshold used to flag unstable LOO diagnostics when
  `criterion = "loo"`. Defaults to `0.7`.

- loo_auto_moment_match:

  Logical; when `TRUE` and `criterion = "loo"`, automatically re-run
  model comparison with `moment_match = TRUE` if any Pareto-k values
  exceed `loo_k_threshold`.

- include_ppc:

  Logical; if `TRUE`, compute \[check_brsm_ppc()\] for both baseline and
  reference models.

- ppc_ndraws:

  Number of posterior predictive draws for PPC.

- ppc_probs:

  Length-2 numeric vector of predictive interval probabilities for PPC.

- seed:

  Optional random seed for reference fitting and PPC.

- data:

  Optional data frame used to fit the reference model when
  `reference_model = NULL` and `object` is not `brsm_fit`.

- response:

  Optional response variable name for automatic reference fitting.
  Inferred from `object` when possible.

- factor_names:

  Optional factor names for automatic reference fitting. Inferred from
  `object` when possible.

- prior:

  Optional prior passed to `brms::brm()` when fitting the reference
  model.

- family:

  Optional family passed to `brms::brm()` when fitting the reference
  model.

- chains:

  Number of chains for automatic reference fitting.

- iter:

  Number of iterations per chain for automatic reference fitting.

- warmup:

  Number of warmup iterations per chain for automatic reference fitting.

- sampling_preset:

  Sampling profile for auto-fitted reference model. One of `"fast"`,
  `"balanced"`, or `"robust"`.

- backend:

  Optional backend passed to `brms::brm()`.

- control:

  Optional named list of NUTS control arguments.

- ...:

  Additional arguments passed to `brms::brm()` when fitting the
  reference model.

## Value

A list with components: `criterion`, `comparison` (from
\[compare_brsm_models()\]), `baseline_model`, `reference_model`,
`reference_type`, `reference_fitted` (logical), optional
`loo_diagnostics` (for `criterion = "loo"`), and optional `ppc` list
with side-by-side summaries.

## Details

The function can either: 1. Accept an already fitted richer model via
`reference_model`, or 2. Fit a richer reference model automatically (for
`brsm_fit` baseline objects, or when fitting inputs are supplied).

## Examples

``` r
if (FALSE) { # \dontrun{
fit_q <- fit_brsm(dat, response = "y", factor_names = c("x1", "x2"))

# Automatic cubic reference fit + LOO comparison
lof <- loftest_brsm(
  object = fit_q,
  reference_type = "cubic",
  include_ppc = TRUE
)
lof$comparison$comparison
lof$ppc$summaries
} # }
```
