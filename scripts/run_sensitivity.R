## ============================================================
## run_sensitivity.R
## Score-quality sensitivity analysis for CALM-VOI.
##
## Holds DGP, n, K fixed and varies the quality of the plug-in
## allocation score sigma_LLM and B_i. Reports VOI plug-in
## variance as a function of score quality, with VOI-oracle as
## the oracle design benchmark.
##
## Sensitivity regimes:
##   O   : oracle sigma (true sigma_LLM, exact)
##   N1  : log-normal noise: sigma_hat = sigma * exp(xi), xi ~ N(0, 0.5^2)
##   N2  : noisier:           xi ~ N(0, 1.0^2)
##   N3  : very noisy:        xi ~ N(0, 1.5^2)
##   R20 : rank-corruption to Spearman rho ~ 0.2
##   R40 : Spearman rho ~ 0.4
##   R60 : Spearman rho ~ 0.6
##   R80 : Spearman rho ~ 0.8
##
## Total cells: 2 DGPs * 8 regimes = 16, all at fixed (n=200, K=5)
## ============================================================

source("calm_voi.R")

## Generate a noisy version of true sigma at given quality regime.
##
## Draws from an isolated RNG stream: the caller's .Random.seed is
## restored on exit, so corrupting the score does not perturb the
## prespecified fold assignment. This lets the eight regimes be run
## as a paired design on identical data and folds.
make_noisy_sigma <- function(sigma_true, regime, rng_state) {
  if (exists(".Random.seed", envir = .GlobalEnv)) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv),
            add = TRUE)
  }
  set.seed(rng_state)
  n <- length(sigma_true)
  if (regime == "O") {
    return(sigma_true)
  } else if (regime == "N1") {
    return(sigma_true * exp(rnorm(n, 0, 0.5)))
  } else if (regime == "N2") {
    return(sigma_true * exp(rnorm(n, 0, 1.0)))
  } else if (regime == "N3") {
    return(sigma_true * exp(rnorm(n, 0, 1.5)))
  } else if (grepl("^R", regime)) {
    ## Rank corruption is only defined when sigma has rank structure.
    ## With constant sigma every quantile is that constant, so the
    ## regime silently returns sigma_true and the cell degenerates
    ## into a re-run of regime O. Fail loudly instead.
    if (sd(sigma_true) < sqrt(.Machine$double.eps)) {
      stop("Regime ", regime, " requires heterogeneous sigma_LLM, but the ",
           "supplied sigma is constant; rank corruption would be a no-op. ",
           "Use a heterogeneous sigma_spec (e.g. \"lognormal\").")
    }
    rho_target <- as.numeric(sub("R", "", regime)) / 100
    ## Blend true ranks with independent random ranks. Weights w and
    ## sqrt(1 - w^2) are orthonormal, so corr(blend, true) = w; the
    ## earlier convex weights (w, 1 - w) gave w/sqrt(w^2 + (1-w)^2),
    ## which overstates the corruption badly (label 0.8 -> rho 0.97).
    w <- rho_target
    ranks_true <- rank(sigma_true)
    ranks_rand <- rank(rnorm(n))
    ranks_blend <- w * ranks_true + sqrt(1 - w^2) * ranks_rand
    ## Map back to a sigma-like scale (preserves the marginal of sigma)
    sigma_hat <- quantile(sigma_true, probs = rank(ranks_blend) / (n + 1),
                          na.rm = TRUE, type = 4, names = FALSE)
    return(sigma_hat)
  } else stop("Unknown regime: ", regime)
}

## Sensitivity DGPs use CONTINUOUS heterogeneous sigma. With the
## constant sigma_spec = 1.0 used previously, sigma_LLM was identical
## across patients, so (a) Neyman degenerated to Unif by construction
## and (b) the four rank-corruption regimes were exact no-ops.
DGPS_SENS <- list(
  "2" = list(dgp = "2", sigma_spec = "lognormal"),
  "3" = list(dgp = "3", sigma_spec = "lognormal")
)
REGIMES <- c("O", "N1", "N2", "N3", "R20", "R40", "R60", "R80")
N_FIXED <- 200
K_FIXED <- 5
B_DEFAULT_SENS <- 5000L
SEED_BASE <- 20260522

## Seed schedule. The grid is a PAIRED design: for a given DGP and
## replication b, every regime sees the same trial data and the same
## cross-fitting fold draws, so cross-regime differences are
## attributable to the allocation score alone rather than to
## independent datasets. Only the corruption draw varies by regime,
## and it is taken from an isolated stream (see make_noisy_sigma).
data_seed_for <- function(dgp_key, b) SEED_BASE + 1e6 * match(dgp_key, names(DGPS_SENS)) + b
fold_seed_for <- function(dgp_key, b) SEED_BASE + 5e6 + 1e6 * match(dgp_key, names(DGPS_SENS)) + b
corr_seed_for <- function(regime, b)  SEED_BASE + 9e6 + 1e4 * match(regime, REGIMES) + b

run_sensitivity_cell <- function(dgp_key, dgp_id, sigma_spec, regime,
                                  n, K, B, seed_offset = NULL) {
  ## We run estimators by hand here because the noised sigma changes
  ## between replications.
  est_names <- c("Unif", "Neyman", "VOI", "VOI_oracle")
  taus <- matrix(NA_real_, nrow = B, ncol = length(est_names),
                 dimnames = list(NULL, est_names))
  ses  <- matrix(NA_real_, nrow = B, ncol = length(est_names),
                 dimnames = list(NULL, est_names))
  rho_realized <- numeric(B)
  t0 <- Sys.time()
  for (b in seq_len(B)) {
    set.seed(data_seed_for(dgp_key, b))
    base   <- generate_baseline(n)
    sigma_eta <- make_sigma_eta(sigma_spec, n)
    out    <- make_outcomes(dgp_id, base$X, base$Tstar, base$A)
    pi_ <- 0.5
    sigma_LLM_true <- 2 * sigma_eta
    B_true <- abs((1 - pi_) * out$beta1 + pi_ * out$beta0) /
              sqrt(pi_ * (1 - pi_))
    q_true <- B_true * sigma_LLM_true

    m_all <- llm_queries(base$X, base$Tstar, sigma_eta, K_max = 25)
    m_first <- m_all[, 1]

    ## Build the noised sigma (isolated RNG stream; leaves the fold
    ## stream below untouched so regimes stay paired)
    sigma_alloc <- make_noisy_sigma(sigma_LLM_true, regime,
                                    rng_state = corr_seed_for(regime, b))
    ## Record the score quality actually realized, so the table can be
    ## labeled by measured rho rather than by nominal target.
    rho_realized[b] <- suppressWarnings(
      cor(sigma_alloc, sigma_LLM_true, method = "spearman"))

    ## One common fold assignment is used by the pilot, allocation,
    ## and every outcome fit. Allocation is performed separately
    ## within held-out folds so one patient's outcome cannot affect
    ## that same patient's budget through a global threshold.
    set.seed(fold_seed_for(dgp_key, b))
    fold_main <- make_folds(n, 5)

    ## Unif
    K_u <- rep(K, n)
    m_bar <- avg_predictions(m_all, K_u)
    r <- aipw_xfit(out$Y, base$A, base$X, M = m_bar, pi_ = pi_,
                   L = 5, include_m = TRUE, fold_id = fold_main)
    taus[b, "Unif"] <- r$tau; ses[b, "Unif"] <- r$se

    ## Neyman (using noised sigma)
    K_n   <- greedy_allocate_by_fold(sigma_alloc, fold_main, K,
                                     K_min = 1L, K_max = 25)
    m_bar <- avg_predictions(m_all, K_n)
    r <- aipw_xfit(out$Y, base$A, base$X, M = m_bar, pi_ = pi_,
                   L = 5, include_m = TRUE,
                   K_alloc = K_n, stratify_by_K = TRUE,
                   fold_id = fold_main)
    taus[b, "Neyman"] <- r$tau; ses[b, "Neyman"] <- r$se

    ## VOI plug-in (B_hat * noised sigma)
    pilot <- cross_fit_pilot_B(base$X, m_first, out$Y, base$A,
                               pi_ = pi_, L = 5, fold_id = fold_main,
                               pilot = "separate")
    q_hat <- pilot$B * sigma_alloc
    K_v   <- greedy_allocate_by_fold(q_hat, fold_main, K,
                                     K_min = 1L, K_max = 25)
    m_bar <- avg_predictions(m_all, K_v)
    r <- aipw_xfit(out$Y, base$A, base$X, M = m_bar, pi_ = pi_,
                   L = 5, include_m = TRUE,
                   K_alloc = K_v, stratify_by_K = TRUE,
                   fold_id = fold_main)
    taus[b, "VOI"] <- r$tau; ses[b, "VOI"] <- r$se

    ## VOI oracle (true B and true sigma)
    K_o   <- greedy_allocate_by_fold(q_true, fold_main, K,
                                     K_min = 1L, K_max = 25)
    m_bar <- avg_predictions(m_all, K_o)
    r <- aipw_xfit(out$Y, base$A, base$X, M = m_bar, pi_ = pi_,
                   L = 5, include_m = TRUE,
                   K_alloc = K_o, stratify_by_K = TRUE,
                   fold_id = fold_main)
    taus[b, "VOI_oracle"] <- r$tau; ses[b, "VOI_oracle"] <- r$se
  }
  dt <- as.numeric(Sys.time() - t0, units = "secs")
  row <- data.frame(DGP = dgp_key, regime = regime, n = n, K = K, B = B,
                    rho_realized = round(mean(rho_realized, na.rm = TRUE), 3),
                    rho_realized_sd = round(sd(rho_realized, na.rm = TRUE), 3),
                    runtime_s = round(dt, 1),
                    stringsAsFactors = FALSE)
  for (e in est_names) {
    bias <- mean(taus[, e]) - TAU_TRUE
    mc_var <- var(taus[, e])
    mc_se_v <- mcse_sample_variance(taus[, e])
    lo <- taus[, e] - 1.96 * ses[, e]
    hi <- taus[, e] + 1.96 * ses[, e]
    cov <- mean(lo <= TAU_TRUE & TAU_TRUE <= hi)
    row[[paste0(e, "_bias")]] <- bias
    row[[paste0(e, "_var")]]  <- mc_var
    row[[paste0(e, "_var_mcse")]] <- mc_se_v
    row[[paste0(e, "_cov")]]  <- cov
  }
  list(row = row, taus = taus, ses = ses)
}

main_sensitivity <- function(B = B_DEFAULT_SENS,
                              save_path = "results_sensitivity.csv",
                              taus_path = "taus_sensitivity.rds",
                              verbose = TRUE, workers = 1L) {
  if (length(B) != 1L || !is.finite(B) || B < 2 || B != round(B)) {
    stop("B must be a finite integer of at least 2")
  }
  B <- as.integer(B)
  if (length(workers) != 1L || !is.finite(workers) ||
      workers < 1 || workers != round(workers)) {
    stop("workers must be a positive integer")
  }
  workers <- as.integer(workers)
  if (.Platform$OS.type == "windows" && workers > 1L) {
    warning("multicore execution is unavailable on Windows; using one worker")
    workers <- 1L
  }
  cells <- expand.grid(dgp_key = names(DGPS_SENS), regime = REGIMES,
                       stringsAsFactors = FALSE)
  run_one <- function(i) {
    dgp_key <- cells$dgp_key[i]
    regime  <- cells$regime[i]
    spec <- DGPS_SENS[[dgp_key]]
    ## Seeds are derived inside run_sensitivity_cell from (DGP, regime, b)
    ## to keep the regimes paired; no per-cell offset is used.
    seed_offset <- NULL
    if (verbose) cat(sprintf("[%d/%d] DGP=%s regime=%s ... ",
                             i, nrow(cells), dgp_key, regime))
    out <- run_sensitivity_cell(dgp_key, spec$dgp, spec$sigma_spec,
                                regime, N_FIXED, K_FIXED, B,
                                seed_offset)
    if (verbose) cat(sprintf("runtime %.0fs\n", out$row$runtime_s))
    list(row = out$row, taus = out$taus, ses = out$ses,
         cell_id = sprintf("%s_%s", dgp_key, regime))
  }
  indices <- seq_len(nrow(cells))
  if (workers > 1L) {
    completed <- parallel::mclapply(
      indices, run_one, mc.cores = workers, mc.preschedule = FALSE
    )
  } else {
    completed <- lapply(indices, run_one)
  }
  rows <- list(); all_taus <- list(); all_ses <- list()
  for (i in indices) {
    out <- completed[[i]]
    if (inherits(out, "try-error")) {
      stop(sprintf("sensitivity cell %d failed: %s", i, out))
    }
    rows[[i]] <- out$row
    all_taus[[out$cell_id]] <- out$taus
    all_ses[[out$cell_id]] <- out$ses
    df <- do.call(rbind, rows)
    write.csv(df, save_path, row.names = FALSE)
    saveRDS(list(taus = all_taus, ses = all_ses), taus_path)
  }
  df <- do.call(rbind, rows)
  write.csv(df, save_path, row.names = FALSE)
  saveRDS(list(taus = all_taus, ses = all_ses), taus_path)
  if (verbose) cat(sprintf("\nDone. Saved %s and %s.\n",
                           save_path, taus_path))
  invisible(df)
}

if (!interactive() && sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  B <- if (length(args) >= 1) as.integer(args[1]) else B_DEFAULT_SENS
  workers <- if (length(args) >= 2) as.integer(args[2]) else 1L
  main_sensitivity(B = B, workers = workers)
}
