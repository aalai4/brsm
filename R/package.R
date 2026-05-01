#' brsm: Bayesian Response Surface Methods
#'
#' @description
#' The \pkg{brsm} package implements Bayesian quadratic response surface models
#' that propagate posterior uncertainty through all optimization quantities.
#' Classical Response Surface Methodology (RSM) provides only point estimates
#' for derived quantities such as stationary points and curvature classification.
#' \pkg{brsm} addresses this by fitting second-order polynomial models via
#' Hamiltonian Monte Carlo (HMC) using \pkg{brms} and \pkg{Stan}, then
#' computing posterior distributions for every optimization-relevant quantity.
#'
#' @section Core Workflow:
#' \enumerate{
#'   \item Fit a Bayesian quadratic model: \code{\link{fit_brsm}}
#'   \item Extract posterior draws: \code{\link{as_brsm_draws}}
#'   \item Compute stationary point: \code{\link{stationary_point}}
#'   \item Classify curvature: \code{\link{canonical_analysis}}
#'   \item Visualize surfaces: \code{\link{brsm_workflow}}
#' }
#'
#' @section Key Functions:
#' \describe{
#'   \item{\code{\link{fit_brsm}}}{Fit a Bayesian quadratic RSM via brms/Stan.}
#'   \item{\code{\link{as_brsm_draws}}}{Extract tidy posterior draws.}
#'   \item{\code{\link{stationary_point}}}{Posterior draws of \eqn{x^* = -\frac{1}{2}B^{-1}b}.}
#'   \item{\code{\link{canonical_analysis}}}{Eigenvalue decomposition of the Hessian \eqn{H = 2B}.}
#'   \item{\code{\link{gradient_quadratic}}}{Posterior gradient \eqn{\nabla f(x) = b + 2Bx}.}
#'   \item{\code{\link{hessian_quadratic}}}{Posterior Hessian \eqn{H = 2B}.}
#'   \item{\code{\link{posterior_ridge_analysis}}}{Bayesian ridge analysis.}
#'   \item{\code{\link{steepest_ascent}}}{Posterior steepest ascent path.}
#'   \item{\code{\link{credible_optimum_region}}}{Credible region for the optimum.}
#'   \item{\code{\link{brsm_workflow}}}{End-to-end analysis and visualization.}
#' }
#'
#' @section References:
#' Alai, A. and Sorwar, A. (2026). \emph{Quantifying Uncertainty in Response
#' Surface Experiments Using Bayesian Quadratic Models.}
#'
#' Bürkner, P.-C. (2017). brms: An R Package for Bayesian Multilevel Models
#' Using Stan. \emph{Journal of Statistical Software}, 80(1), 1--28.
#'
#' Box, G.E.P. and Draper, N.R. (1987). \emph{Empirical Model-Building and
#' Response Surfaces.} Wiley.
#'
#' @import Rcpp
#' @useDynLib brsm, .registration = TRUE
#' @docType package
#' @name brsm-package
#' @aliases brsm
"_PACKAGE"
