# Compute the Posterior Hessian of a Quadratic Response Surface

Returns the Hessian matrix of the quadratic response surface for every
posterior draw. For the quadratic model \\f(x) = \beta_0 + b^\top x +
x^\top B x\\, the Hessian is the constant matrix \\H = 2B\\. Because the
Hessian is constant (does not depend on \\x\\), posterior uncertainty in
\\H\\ arises entirely from posterior uncertainty in the curvature
coefficients \\B\\.

## Usage

``` r
hessian_quadratic(draws, factor_names)
```

## Arguments

- draws:

  A data frame of posterior draws as returned by
  [`as_brsm_draws`](https://aalai4.github.io/brsm/reference/as_brsm_draws.md).

- factor_names:

  Character vector of factor names. Must match columns present in
  `draws`.

## Value

A data frame in long format with one row per (draw, row factor, column
factor) combination containing:

- `draw`:

  Integer index of the posterior draw.

- `row_factor`:

  Factor name for the Hessian row.

- `col_factor`:

  Factor name for the Hessian column.

- `value`:

  The Hessian entry \\H\_{jk} = 2B\_{jk}\\.

## Details

The eigenvalues of \\H\\ determine the curvature class of the stationary
point: all negative \\\Rightarrow\\ maximum, all positive
\\\Rightarrow\\ minimum, mixed \\\Rightarrow\\ saddle point.

## See also

[`gradient_quadratic`](https://aalai4.github.io/brsm/reference/gradient_quadratic.md),
[`canonical_analysis`](https://aalai4.github.io/brsm/reference/canonical_analysis.md),
[`classify_stationary_point`](https://aalai4.github.io/brsm/reference/classify_stationary_point.md)

## Examples

``` r
if (FALSE) { # \dontrun{
draws <- as_brsm_draws(fit)
hess <- hessian_quadratic(draws, factor_names = c("x1", "x2"))
head(hess)

# Posterior mean Hessian matrix
library(dplyr)
hess |>
  group_by(row_factor, col_factor) |>
  summarise(mean_H = mean(value))
} # }
```
