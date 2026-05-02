# Compute the Posterior Gradient of a Quadratic Response Surface

Evaluates the gradient of the fitted quadratic response surface at one
or more points `x` for every posterior draw. The gradient of the
quadratic model \\f(x) = \beta_0 + b^\top x + x^\top B x\\ is \\\nabla
f(x) = b + 2Bx\\, where \\b\\ is the vector of linear coefficients and
\\B\\ is the symmetric curvature matrix.

## Usage

``` r
gradient_quadratic(draws, x, factor_names, normalize = FALSE)
```

## Arguments

- draws:

  A data frame of posterior draws as returned by
  [`as_brsm_draws`](https://aalai4.github.io/brsm/reference/as_brsm_draws.md).

- x:

  A numeric vector (single point) or numeric matrix (multiple points) of
  factor values at which to evaluate the gradient. Columns must match
  `factor_names`.

- factor_names:

  Character vector of factor names. Must match columns present in
  `draws`.

- normalize:

  Logical; if `TRUE`, each gradient vector is scaled to unit length
  (useful for visualizing direction of steepest ascent). Default is
  `FALSE`.

## Value

A data frame with one row per (draw, point) combination containing:

- `draw`:

  Integer index of the posterior draw.

- `point_id`:

  Integer index of the evaluation point.

- factor columns:

  The factor values at which the gradient was evaluated.

- `d_d<factor>`:

  Partial derivative with respect to each factor.

## See also

[`hessian_quadratic`](https://aalai4.github.io/brsm/reference/hessian_quadratic.md),
[`stationary_point`](https://aalai4.github.io/brsm/reference/stationary_point.md),
[`steepest_ascent`](https://aalai4.github.io/brsm/reference/steepest_ascent.md)

## Examples

``` r
if (FALSE) { # \dontrun{
draws <- as_brsm_draws(fit)
# Evaluate gradient at the coded center (0, 0)
grad <- gradient_quadratic(
  draws      = draws,
  x          = c(0, 0),
  factor_names = c("x1", "x2")
)
head(grad)

# Evaluate at multiple points
pts <- matrix(c(-1, 0, 1, -1, 0, 1), ncol = 2)
grad_multi <- gradient_quadratic(draws, x = pts, factor_names = c("x1", "x2"))
} # }
```
