# Decode Stationary Solver Status Codes

Converts per-draw stationary-point solver status codes into an ordered
factor with descriptive labels.

## Usage

``` r
decode_stationary_status(status_code)
```

## Arguments

- status_code:

  Integer-like vector of stationary solver status codes.

## Value

An ordered factor with levels `"ok"`, `"lapack_fail"`,
`"invalid_lu_diag"`, and `"kappa_exceeded"`. Invalid or missing codes
are returned as `NA`.

## Details

Status meanings are:

- `0` = `"ok"`

- `1` = `"lapack_fail"`

- `2` = `"invalid_lu_diag"`

- `3` = `"kappa_exceeded"`
