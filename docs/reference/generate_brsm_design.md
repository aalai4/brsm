# Generate Standard Response Surface Designs

Generates coded experimental designs commonly used in response surface
methodology (RSM): Central Composite Design (CCD) and Box-Behnken Design
(BBD).

## Usage

``` r
generate_brsm_design(
  factor_names,
  design = c("ccd", "bbd"),
  n_center = 4,
  alpha = c("rotatable", "orthogonal"),
  randomize = TRUE,
  seed = NULL,
  output = c("coded", "natural"),
  ranges = NULL
)
```

## Arguments

- factor_names:

  Character vector of factor names.

- design:

  Design type. One of `"ccd"` or `"bbd"`.

- n_center:

  Number of center-point replicates.

- alpha:

  Axial distance used for CCD. Either `"rotatable"`, `"orthogonal"`, or
  a positive numeric scalar.

- randomize:

  Logical; if `TRUE`, randomize run order.

- seed:

  Optional random seed used when `randomize = TRUE`.

- output:

  Output scale. One of `"coded"` (default) or `"natural"`.

- ranges:

  Optional named list of factor ranges in natural units, each a numeric
  vector of length 2. Required when `output = "natural"`.

## Value

A data frame of design runs with columns:

- factor columns

- `.design_type` (ccd/bbd)

- `.block` (factorial/axial/center for CCD; edge/center for BBD)

- `.run_id` (sequential run index)

If `output = "coded"` and `ranges` are supplied, coding metadata is
attached as `brsm_coding` for downstream decode helpers.

## Details

Designs are returned in coded units by default, where factor levels are
centered at 0 and typically span approximately \\\[-1, 1\]\\.
