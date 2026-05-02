# Prepare Factor Columns for BRSM Modeling

Applies lightweight centering/scaling to factor columns for
response-surface modeling and stores coding metadata needed to decode
predictions later.

## Usage

``` r
prepare_brsm_data(
  data,
  factor_names,
  method = c("zscore", "range", "identity")
)
```

## Arguments

- data:

  Data frame containing factor columns.

- factor_names:

  Character vector of factor column names to code.

- method:

  Coding method. One of `"zscore"` (mean/sd) or `"range"`
  (midpoint/half-range; approximately maps observed range to `[-1, 1]`),
  or `"identity"` (no transform; records explicit identity coding
  metadata for already-coded predictors).

## Value

A data frame with coded factor columns. The returned object has class
`brsm_coded_data` and contains attribute `brsm_coding`, a list with
formulas and scaling parameters for each factor.

## Examples

``` r
dat <- data.frame(
  x1 = rnorm(10, 5, 2), x2 = runif(10, 10, 20), y = rnorm(10)
)
coded <- prepare_brsm_data(
  dat, factor_names = c("x1", "x2"), method = "range"
)
str(attr(coded, "brsm_coding"))
#> List of 2
#>  $ method : chr "range"
#>  $ factors:List of 2
#>   ..$ x1:List of 7
#>   .. ..$ center        : num 4.35
#>   .. ..$ scale         : num 3.43
#>   .. ..$ method        : chr "range"
#>   .. ..$ formula       : chr "x1_coded = (x1 - 4.3545831) / 3.4270579"
#>   .. ..$ decode_formula: chr "x1 = x1_coded * 3.4270579 + 4.3545831"
#>   .. ..$ original_range: num [1:2] 0.928 7.782
#>   .. ..$ coded_range   : num [1:2] -1 1
#>   ..$ x2:List of 7
#>   .. ..$ center        : num 14.6
#>   .. ..$ scale         : num 4.07
#>   .. ..$ method        : chr "range"
#>   .. ..$ formula       : chr "x2_coded = (x2 - 14.611968) / 4.0687421"
#>   .. ..$ decode_formula: chr "x2 = x2_coded * 4.0687421 + 14.611968"
#>   .. ..$ original_range: num [1:2] 10.5 18.7
#>   .. ..$ coded_range   : num [1:2] -1 1
```
