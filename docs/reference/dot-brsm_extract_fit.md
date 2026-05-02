# Extract brmsfit from brsm Wrapper

Internal helper that accepts either a `brsm_fit` object or a raw
`brmsfit` and always returns a `brmsfit`.

## Usage

``` r
.brsm_extract_fit(object, caller = "function")
```

## Arguments

- object:

  A `brsm_fit` or `brmsfit` object.

- caller:

  Character scalar used in validation error messages.

## Value

A `brmsfit` object.
