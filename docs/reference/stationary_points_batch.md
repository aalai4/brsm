# Batch Stationary Point Solver

Batch Stationary Point Solver

## Usage

``` r
stationary_points_batch(h_array, b_matrix, kappa_thresh)
```

## Arguments

- h_array:

  Hessian array with dimensions draws x p x p.

- b_matrix:

  Linear coefficient matrix with dimensions draws x p.

- kappa_thresh:

  Condition-number threshold for near-singularity.

## Value

Numeric matrix of stationary points (draws x p), with NA rows for
near-singular or unsolved systems.
