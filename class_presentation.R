# Data Prep ----------------------------------------------------------
library(brsm)
library(brms)   # required for fitting

x=scan(
  text="1 1 1.000 1.000 1.000 6.2 30.5 20.0 0.547 5.06 36.35
 2 1 1.000 -1.000 -1.000 6.2 25.5 5.0 0.643 6.05 21.30
 3 1 -1.000 1.000 -1.000 3.8 30.5 5.0 0.112 1.19 15.96
 4 1 -1.000 -1.000 1.000 3.8 25.5 20.0 0.272 3.96 43.60
 5 1 0.000 0.000 0.000 5.0 28.0 12.5 0.442 3.52 32.74
 6 1 0.000 0.000 0.000 5.0 28.0 12.5 0.108 2.02 19.91
 7 2 1.000 1.000 -1.000 6.2 30.5 5.0 0.752 7.50 79.16
 8 2 1.000 -1.000 1.000 6.2 25.5 20.0 0.615 5.12 23.06
 9 2 -1.000 1.000 1.000 3.8 30.5 20.0 0.063 0.77 18.57
 10 2 -1.000 -1.000 -1.000 3.8 25.5 5.0 0.144 1.25 11.08
 11 2 0.000 0.000 0.000 5.0 28.0 12.5 0.341 3.96 25.29
 12 2 0.000 0.000 0.000 5.0 28.0 12.5 0.315 3.20 46.93
 13 0 -1.633 0.000 0.000 3.0 28.0 12.5 0.637 5.89 59.65
 14 0 1.633 0.000 0.000 7.0 28.0 12.5 1.029 9.98 89.93
 15 0 0.000 -1.633 0.000 5.0 24.0 12.5 0.248 1.85 25.58
 16 0 0.000 1.633 0.000 5.0 32.0 12.5 0.008 0.37 4.21
 17 0 0.000 0.000 -1.633 5.0 28.0 0.0 0.024 0.70 3.72
 18 0 0.000 0.000 1.633 5.0 28.0 25.0 0.037 0.38 8.22
 19 0 0.000 0.000 0.000 5.0 28.0 12.5 0.638 5.74 65.31
 20 0 0.000 0.000 0.000 5.0 28.0 12.5 0.375 3.90 42.27 ",
  what=list(o=0, bl=0, xp=0, xt=0, xg=0, p=0, t=0, gly=0, y=0, spy=0, spa=0),nlines=21)

#obs block xpH xtemp xgly pH temp glycerol yield sp.yld sp.act ;
d=data.frame(x)
attach(d)
names(d)
d$block=factor(bl)
dat <- d[, c("xp", "xt", "xg", "y")]

dat_coded <- prepare_brsm_data(
  data = dat,
  factor_names = c("xp", "xt", "xg"),
  method = "identity"
)

fit <- fit_brsm(
  data = dat_coded,
  response = "y",
  factor_names = c("xp", "xt", "xg"),
  chains = 2,
  iter = 2000,
  warmup = 1000,
  seed = 123,
  sampling_preset = "balanced",
  coding_policy = "ignore", # because already coded
  refresh = 0,
  silent = 2
)

print(fit)
print(summary(fit))
print(check_brsm_fit(fit))
print(check_brsm_ppc(fit))

# Stationary Point Analysis -----------------------------------------------
stat_draws <- stationary_point(object = fit, factor_names = c("xp", "xt", "xg"))

# Summary: posterior mean and credible interval
cat("Posterior mean stationary point:\n")
print(colMeans(stat_draws))

cat("\n95% credible interval (pointwise):\n")
print(apply(stat_draws, 2, quantile, probs = c(0.025, 0.975)))

# Visualization and EDA ---------------------------------------------------
draws_df <- as_brsm_draws(fit, factor_names = c("xp", "xt", "xg")) # Build canonical posterior draws once

p_pairs <- brsm:::plot_posterior_contours(
  draws = draws_df,
  factor_names = c("xp", "xt", "xg"),
  ranges = fit$ranges,
  bins = 12,
  pairwise = TRUE
)
print(p_pairs)