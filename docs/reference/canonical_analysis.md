# Canonical Analysis of Bayesian Response Surface

Computes the eigendecomposition of the quadratic Hessian
(\\\mathbf{B}\\) across posterior draws. Returns principal curvatures
(eigenvalues), canonical axes (eigenvectors), and optional canonical
scores at the stationary point.

## Usage

``` r
canonical_analysis(
  object,
  factor_names = NULL,
  include_scores = TRUE,
  kappa_thresh = 1e+10,
  probs = c(0.025, 0.5, 0.975),
  summary = TRUE
)

# S3 method for class 'brsm_fit'
canonical_analysis(
  object,
  factor_names = NULL,
  include_scores = TRUE,
  kappa_thresh = 1e+10,
  probs = c(0.025, 0.5, 0.975),
  summary = TRUE
)

# S3 method for class 'brmsfit'
canonical_analysis(
  object,
  factor_names = NULL,
  include_scores = TRUE,
  kappa_thresh = 1e+10,
  probs = c(0.025, 0.5, 0.975),
  summary = TRUE
)

# Default S3 method
canonical_analysis(
  object,
  factor_names = NULL,
  include_scores = TRUE,
  kappa_thresh = 1e+10,
  probs = c(0.025, 0.5, 0.975),
  summary = TRUE
)
```

## Arguments

- object:

  A `brsm_fit` object, `brmsfit` object, or data frame of posterior
  draws with Bayesian coefficient columns.

- factor_names:

  Character vector of factor names (required when `object` is a data
  frame).

- include_scores:

  Logical; if `TRUE` (default), compute canonical factor scores at the
  stationary point for each draw.

- kappa_thresh:

  Condition number threshold for Hessian singularity. Draws whose
  Hessian condition number exceeds this are excluded from score
  computation.

- probs:

  Probability levels for posterior summary intervals.

- summary:

  Logical; if `TRUE` (default), return a compact summary data frame with
  posterior means and credible intervals per canonical component. If
  `FALSE`, return per-draw eigenvalues, eigenvectors, and (optionally)
  canonical scores as a list.

## Value

When `summary = TRUE`: a list with components `eigenvalues` (data frame
of posterior summaries for each principal curvature), `eigenvectors`
(data frame of posterior summaries for each axis loading), and
optionally `scores` (data frame of posterior summaries for canonical
factor scores at \\x^\*\\).

When `summary = FALSE`: a list with components `eigenvalues` (draws x p
matrix), `eigenvectors` (draws x p x p array), and optionally `scores`
(draws x p matrix).

## Details

The Hessian \\\mathbf{B}\\ of the quadratic surface \\y = \beta_0 +
\mathbf{b}^\top\mathbf{x} + \mathbf{x}^\top\mathbf{B}\mathbf{x}\\ has
eigendecomposition \\\mathbf{B} =
\mathbf{M}\boldsymbol{\Lambda}\mathbf{M}^\top\\, where the columns of
\\\mathbf{M}\\ are the canonical axes and the diagonal of
\\\boldsymbol{\Lambda}\\ are the principal curvatures.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- fit_brsm(dat, response = "y", factor_names = c("x1", "x2"))
ca <- canonical_analysis(fit)
ca$eigenvalues   # posterior summary of principal curvatures
ca$eigenvectors  # posterior summary of axis loadings
ca$scores        # canonical factor scores at stationary point
} # }
```
