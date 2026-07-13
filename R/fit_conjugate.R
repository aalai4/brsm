#' Exact Conjugate Sampler for Bayesian RSM
#'
#' Fits a second-order polynomial response surface model using the conjugate
#' Normal-Inverse-Gamma prior structure, allowing direct, exact, i.i.d. draws
#' from the joint posterior. Bypasses MCMC sampling and compilation.
#'
#' @keywords internal
.brsm_fit_conjugate <- function(data,
                                response,
                                factor_names,
                                ranges,
                                prior,
                                prior_profile,
                                draws,
                                model_terms) {
  # 1. Build design matrix X and response y
  linear_terms <- factor_names
  interaction_terms <- character(0)
  if (length(factor_names) > 1) {
    combos <- utils::combn(factor_names, 2, simplify = FALSE)
    interaction_terms <- vapply(combos, function(pair) {
      paste0(pair[[1]], ":", pair[[2]])
    }, character(1))
  }
  quadratic_terms <- paste0("I(", factor_names, "^2)")

  rhs_terms <- switch(model_terms,
    first_order = linear_terms,
    first_order_twi = c(linear_terms, interaction_terms),
    pure_quadratic = c(linear_terms, quadratic_terms),
    second_order = c(linear_terms, interaction_terms, quadratic_terms)
  )

  formula_text <- paste(response, "~", paste(rhs_terms, collapse = " + "))
  model_formula <- stats::as.formula(formula_text)

  # Model matrix constructs columns like (Intercept), x1, x2, x1:x2, I(x1^2)
  X <- stats::model.matrix(model_formula, data = data)
  y <- as.numeric(data[[response]])
  n <- length(y)
  p_dims <- ncol(X)

  # 2. Parse/build priors
  if (is.null(prior)) {
    prior <- specify_brsm_priors(
      factor_names = factor_names,
      model_terms = model_terms,
      prior_profile = prior_profile,
      autoscale = TRUE,
      data = data,
      response = response
    )
  }

  prior_df <- as.data.frame(prior)

  # Construct V0 prior covariance diagonal matrix
  V0_diag <- rep(100^2, p_dims) # default wide prior
  names(V0_diag) <- colnames(X)

  # Intercept prior
  intercept_row <- prior_df[prior_df$class == "Intercept", ]
  if (nrow(intercept_row) > 0) {
    V0_diag["(Intercept)"] <- .brsm_parse_prior_sd(intercept_row$prior[1])^2
  }

  # Slope/b priors
  b_global_row <- prior_df[prior_df$class == "b" & (is.na(prior_df$coef) | prior_df$coef == ""), ]
  global_sd <- if (nrow(b_global_row) > 0) .brsm_parse_prior_sd(b_global_row$prior[1]) else 100

  for (col_name in colnames(X)) {
    if (col_name == "(Intercept)") next
    
    # Check for specific coefficient prior
    coef_row <- prior_df[prior_df$class == "b" & prior_df$coef == col_name, ]
    if (nrow(coef_row) > 0) {
      V0_diag[col_name] <- .brsm_parse_prior_sd(coef_row$prior[1])^2
    } else {
      V0_diag[col_name] <- global_sd^2
    }
  }

  V0 <- diag(V0_diag, nrow = p_dims)

  # Sigma prior scale
  sigma_row <- prior_df[prior_df$class == "sigma", ]
  sigma_scale <- if (nrow(sigma_row) > 0) .brsm_parse_prior_sd(sigma_row$prior[1]) else 2.5

  # Inverse-Gamma priors for sigma^2
  a0 <- 1.5
  b0 <- 0.5 * sigma_scale^2

  # 3. Compute Posterior Updates
  V0_inv <- diag(1 / V0_diag, nrow = p_dims)
  XTX <- t(X) %*% X
  V_n_inv <- XTX + V0_inv
  
  # Cholesky decomposition of V_n_inv for numerical stability
  R <- chol(V_n_inv)
  V_n <- chol2inv(R)
  
  XTy <- t(X) %*% y
  m_n <- V_n %*% XTy

  resid <- y - X %*% m_n
  sum_sq_resid <- sum(resid^2)
  prior_penalty <- sum(m_n^2 / V0_diag)
  S_val <- sum_sq_resid + prior_penalty
  
  a_n <- a0 + n / 2
  b_n <- b0 + 0.5 * S_val

  # 4. Generate direct i.i.d draws
  gamma_draws <- stats::rgamma(draws, shape = a_n, rate = b_n)
  sigma2_draws <- 1 / gamma_draws
  sigma_draws <- sqrt(sigma2_draws)

  # draws for beta
  U <- chol(V_n) # U^T U = V_n
  z <- matrix(stats::rnorm(p_dims * draws), nrow = p_dims, ncol = draws)
  scaled_z <- sweep(t(t(U) %*% z), 1, sigma_draws, "*")
  beta_draws <- sweep(scaled_z, 2, m_n, "+")
  colnames(beta_draws) <- colnames(X)

  # standard naming format for draws
  draws_df <- as.data.frame(beta_draws)
  names(draws_df) <- vapply(names(draws_df), function(name) {
    if (name == "(Intercept)") return("b_Intercept")
    if (!startsWith(name, "b_")) return(paste0("b_", name))
    name
  }, character(1))
  draws_df$sigma <- sigma_draws

  # Package output
  fit_obj <- list(
    draws = draws_df,
    data = data,
    X = X,
    y = y,
    V_n = V_n,
    m_n = m_n,
    a_n = a_n,
    b_n = b_n,
    formula = model_formula,
    response = response,
    factor_names = factor_names,
    model_terms = model_terms
  )
  class(fit_obj) <- "brsm_conjugate_fit"
  fit_obj
}

#' Parse Prior Scale/SD from String
#'
#' @keywords internal
.brsm_parse_prior_sd <- function(prior_str) {
  if (is.na(prior_str) || prior_str == "") return(100)
  nums <- as.numeric(regmatches(prior_str, gregexpr("-?[0-9.]+", prior_str))[[1]])
  if (length(nums) == 2) {
    return(nums[2])
  } else if (length(nums) == 3) {
    return(nums[3])
  }
  return(100)
}

#' @export
as.data.frame.brsm_conjugate_fit <- function(x, row.names = NULL, optional = FALSE, ...) {
  x$draws
}

#' @export
summary.brsm_conjugate_fit <- function(object, ...) {
  draws <- object$draws
  coef_cols <- grep("^b_", names(draws), value = TRUE)
  
  fixed_mat <- matrix(NA_real_, nrow = length(coef_cols), ncol = 7)
  rownames(fixed_mat) <- sub("^b_", "", coef_cols)
  colnames(fixed_mat) <- c("Estimate", "Est.Error", "l-95% CI", "u-95% CI", "Rhat", "Bulk_ESS", "Tail_ESS")
  
  for (i in seq_along(coef_cols)) {
    vals <- draws[[coef_cols[i]]]
    fixed_mat[i, "Estimate"] <- mean(vals)
    fixed_mat[i, "Est.Error"] <- stats::sd(vals)
    ci <- stats::quantile(vals, probs = c(0.025, 0.975))
    fixed_mat[i, "l-95% CI"] <- ci[1]
    fixed_mat[i, "u-95% CI"] <- ci[2]
    fixed_mat[i, "Rhat"] <- 1.0
    fixed_mat[i, "Bulk_ESS"] <- nrow(draws)
    fixed_mat[i, "Tail_ESS"] <- nrow(draws)
  }
  
  # Family specific parameters (sigma)
  spec_mat <- matrix(NA_real_, nrow = 1, ncol = 7)
  rownames(spec_mat) <- "sigma"
  colnames(spec_mat) <- c("Estimate", "Est.Error", "l-95% CI", "u-95% CI", "Rhat", "Bulk_ESS", "Tail_ESS")
  
  sig_vals <- draws$sigma
  spec_mat[1, "Estimate"] <- mean(sig_vals)
  spec_mat[1, "Est.Error"] <- stats::sd(sig_vals)
  sig_ci <- stats::quantile(sig_vals, probs = c(0.025, 0.975))
  spec_mat[1, "l-95% CI"] <- sig_ci[1]
  spec_mat[1, "u-95% CI"] <- sig_ci[2]
  spec_mat[1, "Rhat"] <- 1.0
  spec_mat[1, "Bulk_ESS"] <- nrow(draws)
  spec_mat[1, "Tail_ESS"] <- nrow(draws)
  
  out <- list(
    fixed = fixed_mat,
    spec = spec_mat,
    nobs = nrow(object$data)
  )
  class(out) <- "summary_brsm_conjugate"
  out
}

#' @export
print.summary_brsm_conjugate <- function(x, ...) {
  cat("Population-Level Effects:\n")
  print(x$fixed)
  cat("\nFamily Specific Parameters:\n")
  print(x$spec)
}
