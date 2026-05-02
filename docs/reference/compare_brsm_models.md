# Compare Bayesian Models for brsm

Computes model comparison criteria across multiple fitted Bayesian
models.

## Usage

``` r
compare_brsm_models(models, criterion = c("loo", "waic"), ...)
```

## Arguments

- models:

  A named list of `brsm_fit` and/or `brmsfit` objects.

- criterion:

  Comparison criterion: `"loo"` or `"waic"`.

- ...:

  Additional arguments passed to `brms::loo()` or `brms::waic()`.

## Value

A list with components: `criterion`, `estimates` (named list of
criterion objects), and `comparison` (data frame from `loo_compare`).

## Examples

``` r
if (FALSE) { # \dontrun{
fit1 <- fit_brsm(dat, response = "y", factor_names = c("x1", "x2"))
fit2 <- fit_brsm(dat,
  response = "y", factor_names = c("x1", "x2"),
  prior = brms::prior(normal(0, 1), class = "b")
)
cmp <- compare_brsm_models(
  list(default = fit1, stronger_prior = fit2), "loo"
)
cmp$comparison
} # }
```
