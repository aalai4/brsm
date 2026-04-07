% Generated manually; align with R/check_brsm_fit.R
\name{check_brsm_fit}
\alias{check_brsm_fit}
\title{Check MCMC Diagnostics for a brsm Fit}
\usage{
check_brsm_fit(
  object,
  rhat_threshold = 1.01,
  ess_bulk_min = 400,
  ess_tail_min = 400,
  treedepth_limit = NULL,
  bfmi_threshold = 0.3,
  verbose = TRUE
)
}
\arguments{
\item{object}{A \code{brsm_fit} object from \code{fit_brsm()} or a
\code{brmsfit} object.}

\item{rhat_threshold}{Threshold above which Rhat is flagged.}

\item{ess_bulk_min}{Minimum recommended bulk ESS.}

\item{ess_tail_min}{Minimum recommended tail ESS.}

\item{treedepth_limit}{Optional treedepth saturation threshold. If
\code{NULL}, attempts to infer from fit control settings and falls back
to \code{10}.}

\item{bfmi_threshold}{BFMI threshold below which chains are flagged.}

\item{verbose}{Logical; if \code{TRUE}, prints a concise summary.}
}
\value{
A list with components:
\code{overview} (one-row data frame), \code{parameters} (per-parameter
diagnostics), and \code{passed} (logical scalar).
}
\description{
Provides a compact post-fit diagnostics summary focused on fixed-effect
posterior parameters: Rhat, bulk ESS, tail ESS, and NUTS diagnostics
(divergences, treedepth saturation, BFMI).
}
\examples{
\dontrun{
fit <- fit_brsm(dat, response = "y", factor_names = c("x1", "x2"))
diag <- check_brsm_fit(fit)
diag$overview
}
}