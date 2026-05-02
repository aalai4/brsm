# Multi-Response Optimization via Posterior Desirability

Combines response-specific desirability functions over posterior
predictions at candidate points. Predictions are obtained with
\[posterior_predict_brsm()\], mapped to \\\[0, 1\]\\, and aggregated
with a weighted geometric mean.

## Usage

``` r
optimize_brsm_multiresponse(
  models,
  desirability_specs,
  factor_names,
  candidate_points = NULL,
  ranges = NULL,
  n_grid = 25,
  include_residual = FALSE,
  sigma = NULL,
  draw_subset = NULL,
  max_draws = NULL,
  probs = c(0.025, 0.5, 0.975),
  optimize_metric = c("mean", "median"),
  return_draws = FALSE,
  seed = NULL
)
```

## Arguments

- models:

  Named list of response models. Each element can be a `brsm_fit`,
  `brmsfit`, or posterior-draw data frame accepted by
  \[posterior_predict_brsm()\].

- desirability_specs:

  Named list of desirability specifications, one per model. Each element
  must include `goal` and bounds:

  - `goal = "max"` or `"maximize"`: requires `low`, `high`

  - `goal = "min"` or `"minimize"`: requires `low`, `high`

  - `goal = "target"`: requires `low`, `target`, `high`

  Optional fields: `weight` (shape, default 1), `importance`
  (combination weight, default 1).

- factor_names:

  Character vector of factor names.

- candidate_points:

  Optional data frame of candidate points.

- ranges:

  Optional named list of factor ranges used to generate a grid when
  `candidate_points` is `NULL`.

- n_grid:

  Grid resolution per factor when generating candidates from `ranges`.

- include_residual:

  Logical; forwarded to \[posterior_predict_brsm()\].

- sigma:

  Optional residual SD override forwarded to
  \[posterior_predict_brsm()\].

- draw_subset:

  Optional draw subset forwarded to \[posterior_predict_brsm()\].

- max_draws:

  Optional draw cap forwarded to \[posterior_predict_brsm()\].

- probs:

  Probabilities for desirability summaries across draws.

- optimize_metric:

  Which summary metric to maximize when selecting the best point. One of
  `"mean"` or `"median"`.

- return_draws:

  Logical; if `TRUE`, include per-draw combined desirability matrix in
  output.

- seed:

  Optional random seed (used when `include_residual = TRUE`).

## Value

A list with components: `candidate_points` (with desirability
summaries), `best_point` (single-row data frame), `response_at_best`
(per-response summary at the selected point), `draw_count`,
`desirability_specs`, and optional `combined_desirability_draws`.
