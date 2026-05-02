# Check MCMC Diagnostics for a brsm Fit

Provides a compact post-fit diagnostics summary focused on fixed-effect
posterior parameters: Rhat, bulk ESS, tail ESS, and NUTS diagnostics
(divergences, treedepth saturation, BFMI).

## Usage

``` r
check_brsm_fit(
  object,
  rhat_threshold = 1.01,
  ess_bulk_min = 400,
  ess_tail_min = 400,
  treedepth_limit = NULL,
  bfmi_threshold = 0.3,
  verbose = TRUE
)
```

## Arguments

- object:

  A `brsm_fit` object from \[fit_brsm()\] or a `brmsfit` object.

- rhat_threshold:

  Threshold above which Rhat is flagged.

- ess_bulk_min:

  Minimum recommended bulk ESS.

- ess_tail_min:

  Minimum recommended tail ESS.

- treedepth_limit:

  Optional treedepth saturation threshold. If `NULL`, attempts to infer
  from fit control settings and falls back to `10`.

- bfmi_threshold:

  BFMI threshold below which chains are flagged.

- verbose:

  Logical; if `TRUE`, prints a concise summary.

## Value

A list with components: `overview` (one-row data frame), `parameters`
(per-parameter diagnostics), and `passed` (logical scalar).

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- fit_brsm(dat, response = "y", factor_names = c("x1", "x2"))
diag <- check_brsm_fit(fit)
diag$overview
} # }
```
