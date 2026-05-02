# brsm: Bayesian Response Surface Methods

The brsm package implements Bayesian quadratic response surface models
that propagate posterior uncertainty through all optimization
quantities. Classical Response Surface Methodology (RSM) provides only
point estimates for derived quantities such as stationary points and
curvature classification. brsm addresses this by fitting second-order
polynomial models via Hamiltonian Monte Carlo (HMC) using brms and Stan,
then computing posterior distributions for every optimization-relevant
quantity.

## Core Workflow

1.  Fit a Bayesian quadratic model:
    [`fit_brsm`](https://aalai4.github.io/brsm/reference/fit_brsm.md)

2.  Extract posterior draws:
    [`as_brsm_draws`](https://aalai4.github.io/brsm/reference/as_brsm_draws.md)

3.  Compute stationary point:
    [`stationary_point`](https://aalai4.github.io/brsm/reference/stationary_point.md)

4.  Classify curvature:
    [`canonical_analysis`](https://aalai4.github.io/brsm/reference/canonical_analysis.md)

5.  Visualize surfaces:
    [`brsm_workflow`](https://aalai4.github.io/brsm/reference/brsm_workflow.md)

## Key Functions

- [`fit_brsm`](https://aalai4.github.io/brsm/reference/fit_brsm.md):

  Fit a Bayesian quadratic RSM via brms/Stan.

- [`as_brsm_draws`](https://aalai4.github.io/brsm/reference/as_brsm_draws.md):

  Extract tidy posterior draws.

- [`stationary_point`](https://aalai4.github.io/brsm/reference/stationary_point.md):

  Posterior draws of \\x^\* = -\frac{1}{2}B^{-1}b\\.

- [`canonical_analysis`](https://aalai4.github.io/brsm/reference/canonical_analysis.md):

  Eigenvalue decomposition of the Hessian \\H = 2B\\.

- [`gradient_quadratic`](https://aalai4.github.io/brsm/reference/gradient_quadratic.md):

  Posterior gradient \\\nabla f(x) = b + 2Bx\\.

- [`hessian_quadratic`](https://aalai4.github.io/brsm/reference/hessian_quadratic.md):

  Posterior Hessian \\H = 2B\\.

- [`posterior_ridge_analysis`](https://aalai4.github.io/brsm/reference/posterior_ridge_analysis.md):

  Bayesian ridge analysis.

- [`steepest_ascent`](https://aalai4.github.io/brsm/reference/steepest_ascent.md):

  Posterior steepest ascent path.

- [`credible_optimum_region`](https://aalai4.github.io/brsm/reference/credible_optimum_region.md):

  Credible region for the optimum.

- [`brsm_workflow`](https://aalai4.github.io/brsm/reference/brsm_workflow.md):

  End-to-end analysis and visualization.

## References

Alai, A. and Sorwar, A. (2026). *Quantifying Uncertainty in Response
Surface Experiments Using Bayesian Quadratic Models.*

Bürkner, P.-C. (2017). brms: An R Package for Bayesian Multilevel Models
Using Stan. *Journal of Statistical Software*, 80(1), 1–28.

Box, G.E.P. and Draper, N.R. (1987). *Empirical Model-Building and
Response Surfaces.* Wiley.

## See also

Useful links:

- <https://github.com/aalai4/final-project-451-892-a_square>

- Report bugs at
  <https://github.com/aalai4/final-project-451-892-a_square/issues>

## Author

**Maintainer**: Arian Alai <aalai4@nebraska.edu>

Authors:

- Aftab A Sorwar <aftabsorwar15@gmail.com>
