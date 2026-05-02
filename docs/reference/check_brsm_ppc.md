# Posterior Predictive Check Summary for a brsm Fit

Computes a compact posterior predictive summary using simulated
responses from the fitted model.

## Usage

``` r
check_brsm_ppc(
  object,
  ndraws = 200,
  probs = c(0.025, 0.975),
  seed = NULL,
  include_plot = FALSE,
  ...
)
```

## Arguments

- object:

  A `brsm_fit` object from \[fit_brsm()\] or a `brmsfit` object.

- ndraws:

  Number of posterior predictive draws.

- probs:

  Length-2 numeric vector of predictive interval probabilities.

- seed:

  Optional random seed.

- include_plot:

  Logical; if `TRUE`, includes a `pp_check` histogram overlay plot in
  the output.

- ...:

  Additional arguments passed to `brms::posterior_predict()`.

## Value

A list with components: `summary` (one-row data frame), `observed`
(numeric vector), `predicted_mean` (posterior mean prediction per
observation), and optionally `plot`.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- fit_brsm(dat, response = "y", factor_names = c("x1", "x2"))
ppc <- check_brsm_ppc(fit, ndraws = 200)
ppc$summary
} # }
```
