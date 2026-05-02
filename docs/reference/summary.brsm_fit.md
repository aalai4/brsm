# Summary Method for brsm_fit Objects

Display a detailed summary of a fitted Bayesian response surface model,
including MCMC diagnostics and coefficient summaries from the underlying
`brmsfit` object.

## Usage

``` r
# S3 method for class 'brsm_fit'
summary(object, ...)
```

## Arguments

- object:

  An object of class `brsm_fit`.

- ...:

  Additional arguments passed to `brms::summary.brmsfit()`.

## Value

An object of class `summary.brsm_fit` containing: `brsm_fit_obj` (the
original brsm_fit object) and `brmsfit_summary` (summary of the
underlying brmsfit).
