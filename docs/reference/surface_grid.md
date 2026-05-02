# Generate a Regular Grid over Factor Ranges

Creates an expanded prediction grid spanning the specified factor
ranges, which is used internally for surface evaluation and contour
plotting.

## Usage

``` r
surface_grid(ranges, n = 50, center = NULL)
```

## Arguments

- ranges:

  A named list of numeric vectors of length 2, one per factor, giving
  the `c(min, max)` range for each factor.

- n:

  Integer; number of equally-spaced grid points per factor. Must be
  \\\geq 2\\. Default is `50`.

- center:

  Optional named numeric vector. If supplied, this point is guaranteed
  to appear in the grid (appended if not already present).

## Value

A data frame with one column per factor and `n^p` rows (where `p` is the
number of factors), representing the full factorial grid of evaluation
points.
