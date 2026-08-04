## ============================================================
## run_external_pilot.R
## Externally trained frozen-policy experiments.
##
## The rich policy uses a large pilot and a correctly structured basis.
## The reduced policy uses a smaller pilot and a misspecified generic
## interaction basis. The main-effects policy omits leverage
## heterogeneity. Fitted coefficients are frozen before each analysis
## experiment; no analysis outcome or query realization enters its score.
##
## DGPs: 2 (interaction leverage), 3 (arm cancellation)
## This check does not charge the external training cost.
## Analysis n: 200, 400, 800 for the rich policy; 400 otherwise
## Mean query budget: K = 5; cap = 25
## External samples: (n_ext, K_ext) = (5000, 25) or (1000, 5)
##
## Outputs:
##   results_external_pilot.csv
##   results_external_pilot_paired.csv
##   taus_external_pilot.rds
## ============================================================

source("calm_voi.R")

EXTERNAL_ESTIMATORS <- c("Unif", "External", "Oracle")
EXTERNAL_NS <- c(200L, 400L, 800L)
EXTERNAL_DGPS <- c("2", "3")
EXTERNAL_K <- 5L
EXTERNAL_K_MAX <- 25L
EXTERNAL_B_DEFAULT <- 2000L
EXTERNAL_SEED <- 20260711L

EXTERNAL_POLICY_SPECS <- data.frame(
  scenario = c("Rich", "Reduced", "MainOnly"),
  basis = c("rich", "reduced", "main"),
  n_ext = c(5000L, 1000L, 1000L),
  K_ext = c(25L, 5L, 5L),
  stringsAsFactors = FALSE
)

external_design <- function(dgp, X, m, basis = "rich") {
  p <- ncol(X)
  colnames(X) <- paste0("X", seq_len(p))
  if (identical(basis, "main")) {
    return(cbind("(Intercept)" = 1, X, m = m))
  }
  if (identical(basis, "reduced")) {
    xm <- X * m
    colnames(xm) <- paste0("mX", seq_len(p))
    return(cbind("(Intercept)" = 1, X, m = m, xm))
  }
  if (!identical(basis, "rich")) stop("Unknown external-policy basis")
  if (identical(dgp, "2")) {
    xm <- X * m
    colnames(xm) <- paste0("mX", seq_len(p))
    xp <- do.call(cbind, lapply(seq_len(p), function(j) {
      z <- X[, j:p, drop = FALSE] * X[, j]
      colnames(z) <- paste0("X", j, "X", j:p)
      z
    }))
    return(cbind("(Intercept)" = 1, X, m = m, xm, xp))
  }
  if (identical(dgp, "3")) {
    H <- as.numeric(X[, 2] < 0)
    hx <- X * H
    colnames(hx) <- paste0("HX", seq_len(p))
    return(cbind("(Intercept)" = 1, X, H = H, hx, m = m, Hm = H * m))
  }
  stop("The external-pilot experiment is defined only for DGP 2 or 3")
}

fit_external_arm <- function(dgp, X, m, Y, basis) {
  z <- external_design(dgp, X, m, basis)
  theta <- qr.coef(qr(z), Y)
  if (anyNA(theta)) stop("Rank-deficient external calibration fit")
  theta
}

external_derivative <- function(dgp, theta, X, basis) {
  if (identical(basis, "main")) {
    return(rep(unname(theta["m"]), nrow(X)))
  }
  if (identical(basis, "reduced") || identical(dgp, "2")) {
    as.numeric(theta["m"] +
      X %*% theta[paste0("mX", seq_len(ncol(X)))])
  } else {
    H <- as.numeric(X[, 2] < 0)
    as.numeric(theta["m"] + theta["Hm"] * H)
  }
}

fit_frozen_policy <- function(dgp, scenario, basis, n_ext, K_ext) {
  base <- generate_baseline(n_ext)
  out <- make_outcomes(dgp, base$X, base$Tstar, base$A)
  sigma_eta <- rep(1, n_ext)
  queries <- llm_queries(base$X, base$Tstar, sigma_eta, K_ext)
  mbar <- rowMeans(queries)
  theta <- lapply(c(0L, 1L), function(a) {
    ii <- which(base$A == a)
    fit_external_arm(dgp, base$X[ii, , drop = FALSE],
                     mbar[ii], out$Y[ii], basis)
  })
  names(theta) <- c("0", "1")

  beta0_hat <- external_derivative(dgp, theta[["0"]], base$X, basis)
  beta1_hat <- external_derivative(dgp, theta[["1"]], base$X, basis)
  B_hat <- abs(0.5 * beta1_hat + 0.5 * beta0_hat) / 0.5
  B_true <- abs(out$beta1 + out$beta0)
  rho <- suppressWarnings(cor(B_hat, B_true, method = "spearman"))
  if (!is.finite(rho)) rho <- 0
  list(
    dgp = dgp,
    scenario = scenario,
    basis = basis,
    theta = theta,
    n_ext = n_ext,
    K_ext = K_ext,
    calibration_spearman = rho,
    calibration_rmse = sqrt(mean((B_hat - B_true)^2))
  )
}

score_frozen_policy <- function(policy, X, sigma_llm = 2) {
  beta0 <- external_derivative(policy$dgp, policy$theta[["0"]], X,
                               policy$basis)
  beta1 <- external_derivative(policy$dgp, policy$theta[["1"]], X,
                               policy$basis)
  Bhat <- abs(0.5 * beta1 + 0.5 * beta0) / 0.5
  pmax(0, Bhat * sigma_llm)
}

one_external_replication <- function(n, dgp, policy,
                                     K_budget = EXTERNAL_K,
                                     K_max = EXTERNAL_K_MAX) {
  base <- generate_baseline(n)
  out <- make_outcomes(dgp, base$X, base$Tstar, base$A)
  queries <- llm_queries(base$X, base$Tstar, rep(1, n), K_max)
  fold <- make_folds(n, 5)

  q_external <- score_frozen_policy(policy, base$X)
  B_true <- abs(out$beta1 + out$beta0)
  q_oracle <- 2 * B_true
  allocations <- list(
    Unif = rep(K_budget, n),
    External = greedy_allocate(q_external, n * K_budget,
                               K_min = 1L, K_max = K_max),
    Oracle = greedy_allocate(q_oracle, n * K_budget,
                             K_min = 1L, K_max = K_max)
  )

  ans <- lapply(names(allocations), function(method) {
    K <- allocations[[method]]
    mbar <- avg_predictions(queries, K)
    fit <- aipw_xfit(out$Y, base$A, base$X, mbar, L = 5,
                     include_m = TRUE, K_alloc = K,
                     stratify_by_K = method != "Unif", fold_id = fold)
    list(tau = fit$tau, se = fit$se)
  })
  names(ans) <- names(allocations)
  ans
}

summarize_external_cell <- function(taus, ses, dgp, n, B, policy) {
  row <- data.frame(
    scenario = policy$scenario, basis = policy$basis,
    DGP = dgp, n = n, K = EXTERNAL_K, B = B,
    n_ext = policy$n_ext, K_ext = policy$K_ext,
    calibration_spearman = policy$calibration_spearman,
    calibration_rmse = policy$calibration_rmse,
    stringsAsFactors = FALSE
  )
  for (method in EXTERNAL_ESTIMATORS) {
    x <- taus[, method]
    s <- ses[, method]
    row[[paste0(method, "_bias")]] <- mean(x) - TAU_TRUE
    row[[paste0(method, "_var")]] <- var(x)
    row[[paste0(method, "_var_mcse")]] <- mcse_sample_variance(x)
    row[[paste0(method, "_cov")]] <-
      mean(x - 1.96 * s <= TAU_TRUE & TAU_TRUE <= x + 1.96 * s)
    row[[paste0(method, "_mean_se")]] <- mean(s)
  }
  row
}

summarize_external_pair <- function(taus, dgp, n, B, policy) {
  delta <- (taus[, "External"] - TAU_TRUE)^2 -
    (taus[, "Unif"] - TAU_TRUE)^2
  mcse <- sd(delta) / sqrt(B)
  data.frame(
    scenario = policy$scenario, DGP = dgp, n = n, B = B,
    comparison = "External - Unif",
    squared_error_difference = mean(delta), mcse = mcse,
    lower = mean(delta) - 1.96 * mcse,
    upper = mean(delta) + 1.96 * mcse,
    stringsAsFactors = FALSE
  )
}

main_external_pilot <- function(
    B = EXTERNAL_B_DEFAULT,
    save_path = "results_external_pilot.csv",
    paired_path = "results_external_pilot_paired.csv",
    taus_path = "taus_external_pilot.rds",
    verbose = TRUE) {
  if (length(B) != 1L || !is.finite(B) || B < 2 || B != round(B)) {
    stop("B must be a finite integer of at least 2")
  }
  B <- as.integer(B)
  set.seed(EXTERNAL_SEED)
  policies <- list()
  for (j in seq_len(nrow(EXTERNAL_POLICY_SPECS))) {
    spec <- EXTERNAL_POLICY_SPECS[j, ]
    for (dgp in EXTERNAL_DGPS) {
      id <- paste(spec$scenario, dgp, sep = "_")
      policies[[id]] <- fit_frozen_policy(
        dgp = dgp, scenario = spec$scenario, basis = spec$basis,
        n_ext = spec$n_ext, K_ext = spec$K_ext)
    }
  }
  rich_cells <- expand.grid(
    scenario = "Rich", DGP = EXTERNAL_DGPS, n = EXTERNAL_NS,
    stringsAsFactors = FALSE)
  other_cells <- expand.grid(
    scenario = c("Reduced", "MainOnly"), DGP = EXTERNAL_DGPS, n = 400L,
    stringsAsFactors = FALSE)
  cells <- rbind(rich_cells, other_cells)
  rows <- vector("list", nrow(cells))
  paired_rows <- vector("list", nrow(cells))
  raw_taus <- list()
  raw_ses <- list()
  for (i in seq_len(nrow(cells))) {
    dgp <- cells$DGP[i]
    n <- cells$n[i]
    scenario <- cells$scenario[i]
    policy_id <- paste(scenario, dgp, sep = "_")
    taus <- matrix(NA_real_, B, length(EXTERNAL_ESTIMATORS),
                   dimnames = list(NULL, EXTERNAL_ESTIMATORS))
    ses <- taus
    if (verbose) {
      cat(sprintf("[%d/%d] external policy: %s DGP=%s n=%d ... ",
                  i, nrow(cells), scenario, dgp, n))
    }
    t0 <- Sys.time()
    for (b in seq_len(B)) {
      set.seed(EXTERNAL_SEED + 1000000L * i + b)
      z <- one_external_replication(n, dgp, policies[[policy_id]])
      for (method in EXTERNAL_ESTIMATORS) {
        taus[b, method] <- z[[method]]$tau
        ses[b, method] <- z[[method]]$se
      }
    }
    rows[[i]] <- summarize_external_cell(
      taus, ses, dgp, n, B, policies[[policy_id]])
    paired_rows[[i]] <- summarize_external_pair(
      taus, dgp, n, B, policies[[policy_id]])
    id <- sprintf("%s_DGP%s_n%d", scenario, dgp, n)
    raw_taus[[id]] <- taus
    raw_ses[[id]] <- ses
    write.csv(do.call(rbind, rows[seq_len(i)]), save_path,
              row.names = FALSE)
    write.csv(do.call(rbind, paired_rows[seq_len(i)]), paired_path,
              row.names = FALSE)
    saveRDS(list(taus = raw_taus, ses = raw_ses, policies = policies),
            taus_path)
    if (verbose) {
      cat(sprintf("%.1fs\n",
                  as.numeric(Sys.time() - t0, units = "secs")))
    }
  }
  invisible(do.call(rbind, rows))
}

if (!interactive() && sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  B <- if (length(args)) as.integer(args[1]) else EXTERNAL_B_DEFAULT
  main_external_pilot(B = B)
}
