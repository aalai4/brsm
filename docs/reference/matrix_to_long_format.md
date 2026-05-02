# Convert Prediction Matrix to Long-Format Data Frame

Convert Prediction Matrix to Long-Format Data Frame

## Usage

``` r
matrix_to_long_format(pred_matrix, draw_ids, point_ids)
```

## Arguments

- pred_matrix:

  A matrix of predictions (n_draws x n_points)

- draw_ids:

  Integer vector of draw identifiers

- point_ids:

  Integer vector of point identifiers

## Value

A data frame with columns: draw, point_id, estimate
