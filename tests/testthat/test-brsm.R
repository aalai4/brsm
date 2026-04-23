# Clean environment first
rm(list = c("predict_surface", "stationary_point"), inherits = TRUE)

library(brsm)
packageVersion("brsm")

devtools::load_all("c:/Users/Prabhashis-LAPTOP/Downloads/brsm")

get_brsm_fn <- function(name) {
  ns <- asNamespace("brsm")
  if (name %in% getNamespaceExports("brsm")) {
    getExportedValue("brsm", name)
  } else if (exists(name, envir = ns, inherits = FALSE)) {
    get(name, envir = ns, inherits = FALSE)
  } else {
    stop("Function not found in installed brsm: ", name)
  }
}

stationary_point_fn <- get_brsm_fn("stationary_point")
classify_stationarity_point_fn <- get_brsm_fn("classify_stationarity_point")

set.seed(2026)

# Speed/cache settings
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = 2)

# Run controls
N_rep <- 20
DEBUG_METRICS <- FALSE
SAVE_CHECKPOINT <- TRUE

# Clear initial rstan cache
stan_cache_dir <- file.path(tempdir(), "stan_packages")
if (dir.exists(stan_cache_dir)) {
  try({
    dll_files <- list.files(
      stan_cache_dir,
      pattern = "\\.dll$",
      full.names = TRUE
    )
    for (dll in dll_files) {
      tryCatch(dyn.unload(dll), error = function(e) NULL)
    }
  }, silent = TRUE)
}

# -----------------------------------------------------------------------------
# Helpers: truth
# -----------------------------------------------------------------------------
make_true_params <- function(p) {
  if (p == 2) {
    beta0 <- 5
    b <- c(2, 1)
    B <- matrix(c(
      -1, 0.2,
      0.2, -0.8
    ), nrow = 2, byrow = TRUE)
  } else {
    beta0 <- 5
    b <- c(2, 1, -0.5)
    B <- matrix(c(
      -1,   0.2,  0.1,
      0.2, -0.8, 0.15,
      0.1,  0.15, -1.2
    ), nrow = 3, byrow = TRUE)
  }
  list(beta0 = beta0, b = b, B = B)
}

true_stationary_point <- function(params) {
  -0.5 * solve(params$B) %*% params$b
}

true_curvature_class <- function(params) {
  eigs <- eigen(2 * params$B, symmetric = TRUE, only.values = TRUE)$values
  if (all(eigs < 0)) "maximum" else if (all(eigs > 0)) "minimum" else "saddle"
}

# -----------------------------------------------------------------------------
# Helper: synthetic dataset
# -----------------------------------------------------------------------------
make_dataset <- function(n, p, params, sigma) {
  factor_names <- paste0("x", seq_len(p))
  design <- as.data.frame(
    matrix(
      runif(n * p, -1, 1),
      nrow = n, ncol = p,
      dimnames = list(NULL, factor_names)
    )
  )
  X <- as.matrix(design)
  eta <- params$beta0 + X %*% params$b + rowSums((X %*% params$B) * X)
  y <- as.numeric(eta) + rnorm(n, 0, sigma)
  data.frame(y = y, design)
}

# -----------------------------------------------------------------------------
# Helper: robust draw extraction (handles names like b_Ix1E2)
# -----------------------------------------------------------------------------
extract_draws_robust <- function(fit_obj, factor_names, debug = FALSE) {
  raw <- as.data.frame(fit_obj$fit, check.names = FALSE)
  raw_cols <- names(raw)

  find_first <- function(cands) {
    h <- intersect(cands, raw_cols)
    if (length(h) == 0) NA_character_ else h[[1]]
  }

  out <- list()

  c_int <- find_first(c("b_Intercept", "b_Intercept."))
  if (is.na(c_int)) stop("Could not find intercept column.")
  out[["b_Intercept"]] <- raw[[c_int]]

  for (f in factor_names) {
    c_lin <- find_first(c(paste0("b_", f), paste0("b_", f, ".")))
    if (is.na(c_lin)) stop("Could not find linear column for ", f)
    out[[paste0("b_", f)]] <- raw[[c_lin]]

    c_quad <- find_first(c(
      paste0("b_I(", f, "^2)"),
      paste0("b_I", f, "E2"),
      paste0("b_I.", f, ".2."),
      paste0("b_I.", f, ".2"),
      paste0("b_I", f, ".2."),
      paste0("b_I", f, ".2"),
      make.names(paste0("b_I(", f, "^2)"))
    ))

    if (is.na(c_quad)) {
      idx <- grep(
        paste0("^b_?I.*", f, ".*(\\^2|E2|\\.2\\.?|2\\)?)$"),
        raw_cols,
        perl = TRUE
      )
      if (length(idx) > 0) c_quad <- raw_cols[idx[1]]
    }

    if (is.na(c_quad)) {
      stop(
        "Could not find quadratic column for ", f, ". Raw b_ names: ",
        paste(grep("^b_", raw_cols, value = TRUE), collapse = ", ")
      )
    }

    out[[paste0("b_I(", f, "^2)")]] <- raw[[c_quad]]
  }

  if (length(factor_names) > 1) {
    for (i in seq_len(length(factor_names) - 1)) {
      for (j in (i + 1):length(factor_names)) {
        f1 <- factor_names[i]
        f2 <- factor_names[j]
        c_intx <- find_first(c(
          paste0("b_", f1, ":", f2),
          paste0("b_", f2, ":", f1),
          paste0("b_", f1, ".", f2),
          paste0("b_", f2, ".", f1)
        ))
        if (!is.na(c_intx)) out[[paste0("b_", f1, ":", f2)]] <- raw[[c_intx]]
      }
    }
  }

  draws <- as.data.frame(out, check.names = FALSE)

  if (debug) {
    cat("raw b_ names:\n")
    print(grep("^b_", raw_cols, value = TRUE))
    cat("final mapped names:\n")
    print(names(draws))
  }

  draws
}

safe_mean <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

# -----------------------------------------------------------------------------
# Helper: one replicate
# -----------------------------------------------------------------------------
run_replicate <- function(
  data,
  factor_names,
  true_x_star,
  true_class,
  beta_true_vec,
  debug_metrics = FALSE) {
  fit <- tryCatch(
    fit_brsm(
      data = data,
      response = "y",
      factor_names = factor_names,
      chains = 2,
      iter = 1000,
      warmup = 500,
      seed = sample.int(.Machine$integer.max, 1),
      sampling_preset = "fast",
      refresh = 0,
      silent = 2,
      cores = 2
    ),
    error = function(e) {
      message("fit error: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(fit)) return(NULL)

  draws <- tryCatch(
    extract_draws_robust(fit, factor_names, debug = FALSE),
    error = function(e) {
      message("draw extraction error: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(draws)) return(NULL)

  coef_names <- c(
    "b_Intercept",
    paste0("b_", factor_names),
    paste0("b_I(", factor_names, "^2)")
  )
  if (!all(coef_names %in% names(draws))) {
    message(
      "missing coef columns: ",
      paste(setdiff(coef_names, names(draws)), collapse = ", ")
    )
    return(NULL)
  }

  coef_draws <- draws[, coef_names, drop = FALSE]

  post_means <- colMeans(coef_draws)
  bias <- mean(post_means - beta_true_vec)
  rmse <- sqrt(mean((post_means - beta_true_vec)^2))

  ci_lo <- apply(coef_draws, 2, quantile, probs = 0.025, na.rm = TRUE)
  ci_hi <- apply(coef_draws, 2, quantile, probs = 0.975, na.rm = TRUE)
  coverage <- mean(beta_true_vec >= ci_lo & beta_true_vec <= ci_hi)

  # Stationary point
  sp <- tryCatch(
    stationary_point_fn(draws, factor_names = factor_names, kappa_thresh = Inf),
    error = function(e) {
      if (debug_metrics) {
        message("stationary_point error: ", conditionMessage(e))
      }
      NULL
    }
  )

  sp_error <- if (!is.null(sp)) {
    sp_ok <- sp[stats::complete.cases(sp), , drop = FALSE]
    if (nrow(sp_ok) > 0) {
      sqrt(sum((colMeans(sp_ok) - as.numeric(true_x_star))^2))
    } else {
      NA_real_
    }
  } else {
    NA_real_
  }

  # Classification
  cl <- tryCatch(
    classify_stationarity_point_fn(draws, factor_names = factor_names),
    error = function(e) {
      if (debug_metrics) {
        message("classify_stationarity_point error: ", conditionMessage(e))
      }
      NULL
    }
  )

  class_prob <- if (!is.null(cl) && "classification" %in% names(cl)) {
    cls <- as.character(cl$classification)
    valid_cls <- !is.na(cls)
    if (any(valid_cls)) mean(cls[valid_cls] == true_class) else NA_real_
  } else {
    NA_real_
  }

  if (debug_metrics) {
    sp_ok_n <- if (!is.null(sp)) sum(stats::complete.cases(sp)) else 0
    cat("    debug: sp complete rows =", sp_ok_n, "\n")
    if (!is.null(cl) && "classification" %in% names(cl)) {
      print(table(as.character(cl$classification), useNA = "ifany"))
    } else {
      cat("    debug: classification unavailable\n")
    }
  }

  list(
    bias = bias,
    rmse = rmse,
    coverage = coverage,
    sp_error = sp_error,
    class_prob = class_prob
  )
}

# -----------------------------------------------------------------------------
# Simulation grid
# -----------------------------------------------------------------------------
configs <- expand.grid(
  p = c(2, 3),
  n = c(20, 50, 100),
  sigma = c(0.5, 1.0, 2.0),
  stringsAsFactors = FALSE
)
configs <- configs[order(configs$p, configs$n, configs$sigma), , drop = FALSE]

results <- vector("list", nrow(configs))

for (cfg_i in seq_len(nrow(configs))) {
  p <- configs$p[cfg_i]
  n <- configs$n[cfg_i]
  sigma <- configs$sigma[cfg_i]

  factor_names <- paste0("x", seq_len(p))
  params <- make_true_params(p)
  true_x_star <- true_stationary_point(params)
  true_class <- true_curvature_class(params)
  beta_true_vec <- c(params$beta0, params$b, diag(params$B))

  cat(sprintf(
    "\nConfig %d/%d: p=%d, n=%d, sigma=%.1f [true class: %s]\n",
    cfg_i, nrow(configs), p, n, sigma, true_class
  ))

  # ==========================================================================
  # PERIODIC CACHE CLEANUP (before each config after the first)
  # ==========================================================================
  if (cfg_i > 1) {
    cat("  [Clearing rstan DLL cache...]\n")
    gc()  # Force garbage collection

    # Unload cached stan DLLs from temp directory
    stan_cache_dir <- file.path(tempdir(), "stan_packages")
    if (dir.exists(stan_cache_dir)) {
      tryCatch({
        dll_files <- list.files(
          stan_cache_dir,
          pattern = "\\.dll$",
          full.names = TRUE
        )
        for (dll in dll_files) {
          tryCatch(dyn.unload(dll), error = function(e) NULL)
        }
      }, error = function(e) NULL)
    }
    Sys.sleep(0.5)  # Give Windows time to release locks
  }
  # ==========================================================================

  rep_results <- vector("list", N_rep)

  for (r in seq_len(N_rep)) {
    if (r %% 10 == 0) cat(sprintf("  replicate %d/%d\n", r, N_rep))
    data_r <- make_dataset(n, p, params, sigma)
    rep_results[[r]] <- run_replicate(
      data = data_r,
      factor_names = factor_names,
      true_x_star = true_x_star,
      true_class = true_class,
      beta_true_vec = beta_true_vec,
      debug_metrics = DEBUG_METRICS
    )
  }

  valid <- Filter(Negate(is.null), rep_results)
  n_valid <- length(valid)

  cat(sprintf("  n_valid = %d/%d\n", n_valid, N_rep))

  if (n_valid > 0) {
    results[[cfg_i]] <- data.frame(
      p = p,
      n = n,
      sigma = sigma,
      n_valid = n_valid,
      bias = safe_mean(sapply(valid, "[[", "bias")),
      rmse = safe_mean(sapply(valid, "[[", "rmse")),
      coverage = safe_mean(sapply(valid, "[[", "coverage")),
      sp_error = safe_mean(sapply(valid, "[[", "sp_error")),
      class_prob = safe_mean(sapply(valid, "[[", "class_prob"))
    )
  } else {
    results[[cfg_i]] <- data.frame(
      p = p, n = n, sigma = sigma, n_valid = 0,
      bias = NA_real_, rmse = NA_real_, coverage = NA_real_,
      sp_error = NA_real_, class_prob = NA_real_
    )
  }

  cat(sprintf(
    "  Done: coverage=%.3f, sp_error=%.3f, class_prob=%.3f\n",
    results[[cfg_i]]$coverage,
    results[[cfg_i]]$sp_error,
    results[[cfg_i]]$class_prob
  ))

  if (SAVE_CHECKPOINT) {
    saveRDS(do.call(rbind, results[seq_len(cfg_i)]), "sim_results_partial.rds")
  }
}

# Compile final results
results_final <- do.call(rbind, Filter(function(x) !is.null(x), results))
rownames(results_final) <- NULL

cat("\n\n===== FINAL RESULTS =====\n")
print(results_final)

saveRDS(results_final, "sim_results_final.rds")
