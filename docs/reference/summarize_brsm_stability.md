# Summarize BRSM Stability Diagnostics

Extracts and displays stability diagnostics from BRSM analysis results.
Provides a unified view of exclusion rates, failure modes, and numerical
stability across function calls.

## Usage

``` r
summarize_brsm_stability(object, ...)

# S3 method for class 'brsm_fit'
summarize_brsm_stability(object, ...)

# Default S3 method
summarize_brsm_stability(object, ...)
```

## Arguments

- object:

  A \`brsm_fit\` object, or a result object (data frame or list) from
  \`stationary_point()\`, \`posterior_ridge_analysis()\`,
  \`canonical_analysis()\`, or similar functions that attach a
  \`"diagnostics"\` attribute.

- ...:

  Additional result objects to summarize together.

## Value

A data frame with one row per analyzed result, containing:

- \`function_name\`: Name of the analysis function.

- \`n_draws\`: Total number of posterior draws.

- \`n_excluded\`: Number of excluded draws (failed status).

- \`pct_excluded\`: Percentage of excluded draws.

- \`status_ok\`: Count of OK solves.

- \`status_lapack_fail\`: Count of LAPACK failures.

- \`status_invalid_lu_diag\`: Count of invalid LU diagonal failures.

- \`status_kappa_exceeded\`: Count of high-condition-number exclusions.

## Details

For \`brsm_fit\` objects, this function extracts diagnostics from all
major internal analyses (stationary point computation, etc.).

For analysis result objects with a \`"diagnostics"\` attribute, it
extracts the status codes and aggregates them into counts.

Use this to quickly assess whether numerical instability is affecting
your inferences, and to understand which types of failures are most
prevalent.
