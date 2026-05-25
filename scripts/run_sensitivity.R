## ============================================================
## run_sensitivity.R
## Score-quality sensitivity analysis for CALM-VOI.
##
## Holds DGP, n, K fixed and varies the quality of the plug-in
## allocation score sigma_LLM and B_i. Reports VOI plug-in
## variance as a function of score quality, with VOI-oracle as
## upper-bound benchmark.
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

## Robust path resolution: works from repo root or scripts/
.calm_voi_path <- if (file.exists("R/calm_voi.R")) "R/calm_voi.R" else
                  if (file.exists("../R/calm_voi.R")) "../R/calm_voi.R" else
                  stop("Cannot find R/calm_voi.R. Run from repo root or scripts/.")
source(.calm_voi_path)


if (!dir.exists("outputs")) dir.create("outputs", recursive = TRUE)

## Generate a noisy version of true sigma at given quality regime
make_noisy_sigma <- function(sigma_true, regime, rng_state) {
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
    rho_target <- as.numeric(sub("R", "", regime)) / 100
    ## Generate a noisy version with approximate Spearman rho
    ## by blending true ranks with random ranks
    ranks_true <- rank(sigma_true)
    ranks_rand <- rank(rnorm(n))
    ## Convex combination of rank vectors (approximate, not exact rho)
    w <- rho_target
    ranks_blend <- w * ranks_true + (1 - w) * ranks_rand
    ## Map back to a sigma-like scale (median = median of true sigma)
    sigma_hat <- quantile(sigma_true, probs = rank(ranks_blend) / (n + 1),
                          na.rm = TRUE, type = 4, names = FALSE)
    return(sigma_hat)
  } else stop("Unknown regime: ", regime)
}

DGPS_SENS <- list(
  "2" = list(dgp = "2", sigma_spec = 1.0),
  "3" = list(dgp = "3", sigma_spec = 1.0)
)
REGIMES <- c("O", "N1", "N2", "N3", "R20", "R40", "R60", "R80")
N_FIXED <- 200
K_FIXED <- 5
B_DEFAULT_SENS <- 1000
SEED_BASE <- 20260522

run_sensitivity_cell <- function(dgp_key, dgp_id, sigma_spec, regime,
                                  n, K, B, seed_offset) {
  ## We run estimators by hand here because the noised sigma changes
  ## between replications.
  est_names <- c("Unif", "Neyman", "VOI", "VOI_oracle")
  taus <- matrix(NA_real_, nrow = B, ncol = length(est_names),
                 dimnames = list(NULL, est_names))
  ses  <- matrix(NA_real_, nrow = B, ncol = length(est_names),
                 dimnames = list(NULL, est_names))
  t0 <- Sys.time()
  for (b in seq_len(B)) {
    set.seed(seed_offset + b)
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

    ## Build the noised sigma
    sigma_alloc <- make_noisy_sigma(sigma_LLM_true, regime,
                                    rng_state = seed_offset + b + 1e6)

    ## Unif
    K_u <- rep(K, n)
    m_bar <- avg_predictions(m_all, K_u)
    r <- aipw_xfit(out$Y, base$A, base$X, M = m_bar, pi_ = pi_,
                   L = 5, include_m = TRUE)
    taus[b, "Unif"] <- r$tau; ses[b, "Unif"] <- r$se

    ## Neyman (using noised sigma)
    K_n_f <- water_filling(sigma_alloc, n * K, K_min = 1)
    K_n   <- integer_allocate(K_n_f, n * K, K_min = 1, K_max_cap = 25)
    m_bar <- avg_predictions(m_all, K_n)
    r <- aipw_xfit(out$Y, base$A, base$X, M = m_bar, pi_ = pi_,
                   L = 5, include_m = TRUE,
                   K_alloc = K_n, stratify_by_K = TRUE)
    taus[b, "Neyman"] <- r$tau; ses[b, "Neyman"] <- r$se

    ## VOI plug-in (B_hat * noised sigma)
    pilot <- cross_fit_pilot_B(base$X, m_first, out$Y, base$A,
                               pi_ = pi_, L = 5)
    q_hat <- pilot$B * sigma_alloc
    K_v_f <- water_filling(q_hat, n * K, K_min = 1)
    K_v   <- integer_allocate(K_v_f, n * K, K_min = 1, K_max_cap = 25)
    m_bar <- avg_predictions(m_all, K_v)
    r <- aipw_xfit(out$Y, base$A, base$X, M = m_bar, pi_ = pi_,
                   L = 5, include_m = TRUE,
                   K_alloc = K_v, stratify_by_K = TRUE)
    taus[b, "VOI"] <- r$tau; ses[b, "VOI"] <- r$se

    ## VOI oracle (true B and true sigma)
    K_o_f <- water_filling(q_true, n * K, K_min = 1)
    K_o   <- integer_allocate(K_o_f, n * K, K_min = 1, K_max_cap = 25)
    m_bar <- avg_predictions(m_all, K_o)
    r <- aipw_xfit(out$Y, base$A, base$X, M = m_bar, pi_ = pi_,
                   L = 5, include_m = TRUE,
                   K_alloc = K_o, stratify_by_K = TRUE)
    taus[b, "VOI_oracle"] <- r$tau; ses[b, "VOI_oracle"] <- r$se
  }
  dt <- as.numeric(Sys.time() - t0, units = "secs")
  row <- data.frame(DGP = dgp_key, regime = regime, n = n, K = K, B = B,
                    runtime_s = round(dt, 1),
                    stringsAsFactors = FALSE)
  for (e in est_names) {
    bias <- mean(taus[, e]) - TAU_TRUE
    mc_var <- var(taus[, e])
    mc_se_v <- sqrt(mc_var * (2 / (B - 1)))
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
                              save_path = "outputs/results_sensitivity.csv",
                              taus_path = "outputs/taus_sensitivity.rds",
                              verbose = TRUE) {
  cells <- expand.grid(dgp_key = names(DGPS_SENS), regime = REGIMES,
                       stringsAsFactors = FALSE)
  rows <- list(); all_taus <- list(); all_ses <- list()
  for (i in seq_len(nrow(cells))) {
    dgp_key <- cells$dgp_key[i]
    regime  <- cells$regime[i]
    spec <- DGPS_SENS[[dgp_key]]
    seed_offset <- SEED_BASE + 1e6 * (1000 + i)
    if (verbose) cat(sprintf("[%d/%d] DGP=%s regime=%s ... ",
                             i, nrow(cells), dgp_key, regime))
    out <- run_sensitivity_cell(dgp_key, spec$dgp, spec$sigma_spec,
                                regime, N_FIXED, K_FIXED, B,
                                seed_offset)
    rows[[i]] <- out$row
    all_taus[[sprintf("%s_%s", dgp_key, regime)]] <- out$taus
    all_ses[[sprintf("%s_%s", dgp_key, regime)]]  <- out$ses
    if (verbose) cat(sprintf("runtime %.0fs\n", out$row$runtime_s))
    df <- do.call(rbind, rows)
    write.csv(df, save_path, row.names = FALSE)
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
  main_sensitivity(B = B)
}
