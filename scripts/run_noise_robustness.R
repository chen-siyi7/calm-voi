## ============================================================
## run_noise_robustness.R
## Query-noise stress test outside the Gaussian-iid working model.
##
## DGP 3, n = 400, mean K = 5, cap = 25.
## Regimes:
##   Gaussian       independent standard Gaussian query errors
##   Skewed         independent centered, standardized lognormal errors
##   Correlated     Gaussian errors with within-patient correlation 0.30
##
## Outputs:
##   results_noise_robustness.csv
##   taus_noise_robustness.rds
## ============================================================

source("calm_voi.R")

NOISE_REGIMES <- c("Gaussian", "Skewed", "Correlated")
NOISE_ESTIMATORS <- c("Unif", "VOI_oracle")
NOISE_N <- 400L
NOISE_K <- 5L
NOISE_K_MAX <- 25L
NOISE_RHO <- 0.30
NOISE_B_DEFAULT <- 2000L
NOISE_SEED <- 20260719L

standardized_lognormal <- function(n, sdlog = 0.75) {
  mean_x <- exp(sdlog^2 / 2)
  sd_x <- sqrt((exp(sdlog^2) - 1) * exp(sdlog^2))
  (rlnorm(n, meanlog = 0, sdlog = sdlog) - mean_x) / sd_x
}

robust_queries <- function(X, Tstar, regime, K_max = NOISE_K_MAX,
                           rho = NOISE_RHO) {
  n <- nrow(X)
  if (identical(regime, "Gaussian")) {
    eta <- matrix(rnorm(n * K_max), nrow = n)
  } else if (identical(regime, "Skewed")) {
    eta <- matrix(standardized_lognormal(n * K_max), nrow = n)
  } else if (identical(regime, "Correlated")) {
    common <- matrix(rnorm(n), nrow = n, ncol = K_max)
    idio <- matrix(rnorm(n * K_max), nrow = n)
    eta <- sqrt(rho) * common + sqrt(1 - rho) * idio
  } else {
    stop("Unknown query-noise regime: ", regime)
  }
  xb <- as.numeric(X %*% BETA_X)
  2 * (matrix(Tstar, n, K_max) + eta) + matrix(xb, n, K_max)
}

one_noise_replication <- function(n = NOISE_N, K_budget = NOISE_K) {
  base <- generate_baseline(n)
  out <- make_outcomes("3", base$X, base$Tstar, base$A)
  fold <- make_folds(n, 5)
  q_oracle <- 2 * abs(out$beta1 + out$beta0)
  K_oracle <- greedy_allocate(q_oracle, n * K_budget,
                              K_min = 1L, K_max = NOISE_K_MAX)
  allocations <- list(
    Unif = rep(K_budget, n),
    VOI_oracle = K_oracle
  )
  ans <- list()
  for (regime in NOISE_REGIMES) {
    queries <- robust_queries(base$X, base$Tstar, regime)
    ans[[regime]] <- lapply(names(allocations), function(method) {
      K <- allocations[[method]]
      mbar <- avg_predictions(queries, K)
      fit <- aipw_xfit(out$Y, base$A, base$X, mbar, L = 5,
                       include_m = TRUE, K_alloc = K,
                       stratify_by_K = method != "Unif",
                       fold_id = fold)
      list(tau = fit$tau, se = fit$se)
    })
    names(ans[[regime]]) <- names(allocations)
  }
  ans
}

main_noise_robustness <- function(
    B = NOISE_B_DEFAULT,
    save_path = "results_noise_robustness.csv",
    taus_path = "taus_noise_robustness.rds",
    verbose = TRUE) {
  if (length(B) != 1L || !is.finite(B) || B < 2 || B != round(B)) {
    stop("B must be a finite integer of at least 2")
  }
  B <- as.integer(B)
  taus <- lapply(NOISE_REGIMES, function(z)
    matrix(NA_real_, B, length(NOISE_ESTIMATORS),
           dimnames = list(NULL, NOISE_ESTIMATORS)))
  ses <- taus
  names(taus) <- names(ses) <- NOISE_REGIMES
  t0 <- Sys.time()
  for (b in seq_len(B)) {
    set.seed(NOISE_SEED + b)
    z <- one_noise_replication()
    for (regime in NOISE_REGIMES) {
      for (method in NOISE_ESTIMATORS) {
        taus[[regime]][b, method] <- z[[regime]][[method]]$tau
        ses[[regime]][b, method] <- z[[regime]][[method]]$se
      }
    }
    if (verbose && b %% max(1L, B %/% 10L) == 0L) {
      cat(sprintf("noise robustness: %d/%d\n", b, B))
    }
  }
  rows <- lapply(NOISE_REGIMES, function(regime) {
    row <- data.frame(
      DGP = "3", n = NOISE_N, K = NOISE_K, B = B,
      regime = regime,
      rho = if (regime == "Correlated") NOISE_RHO else 0,
      stringsAsFactors = FALSE
    )
    for (method in NOISE_ESTIMATORS) {
      x <- taus[[regime]][, method]
      s <- ses[[regime]][, method]
      row[[paste0(method, "_bias")]] <- mean(x) - TAU_TRUE
      row[[paste0(method, "_var")]] <- var(x)
      row[[paste0(method, "_var_mcse")]] <- mcse_sample_variance(x)
      row[[paste0(method, "_cov")]] <-
        mean(x - 1.96 * s <= TAU_TRUE & TAU_TRUE <= x + 1.96 * s)
      row[[paste0(method, "_mean_se")]] <- mean(s)
    }
    row
  })
  df <- do.call(rbind, rows)
  write.csv(df, save_path, row.names = FALSE)
  saveRDS(list(taus = taus, ses = ses), taus_path)
  if (verbose) {
    cat(sprintf("Done in %.1fs. Saved %s and %s.\n",
                as.numeric(Sys.time() - t0, units = "secs"),
                save_path, taus_path))
  }
  invisible(df)
}

if (!interactive() && sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  B <- if (length(args)) as.integer(args[1]) else NOISE_B_DEFAULT
  main_noise_robustness(B = B)
}
