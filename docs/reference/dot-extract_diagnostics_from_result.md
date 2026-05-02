# Extract Diagnostics from a Result Object

Internal helper to extract stability diagnostics from an analysis
result.

## Usage

``` r
.extract_diagnostics_from_result(obj)
```

## Arguments

- obj:

  A result object with a \`"diagnostics"\` attribute.

## Value

A one-row data frame with aggregated diagnostics, or NULL if no
diagnostics are found.
