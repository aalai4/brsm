# Detailed Batch Stationary Point Solver

Detailed Batch Stationary Point Solver

## Usage

``` r
stationary_points_batch_details(h_array, b_matrix, kappa_thresh)
```

## Arguments

- h_array:

  Hessian array with dimensions draws x p x p.

- b_matrix:

  Linear coefficient matrix with dimensions draws x p.

- kappa_thresh:

  Condition-number threshold for near-singularity.

## Value

A list with stationary points, per-draw status codes, and kappa proxies.
