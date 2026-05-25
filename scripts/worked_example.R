## ============================================================
## worked_example.R
## A semi-synthetic worked example illustrating CALM-VOI on a
## simulated symptomatic-treatment trial.
##
## IMPORTANT FRAMING:
##   This is NOT a real-data application. It is fully synthetic
##   data designed to mimic the structure of a 250-patient
##   clinical trial. The "LLM" is a simple multilayer perceptron
##   pre-trained on synthetic training data and queried with
##   per-query Gaussian noise (mimicking temperature sampling).
##
##   A genuine real-data application requires:
##     - Actual patient text records (e.g., MIMIC-III clinical notes)
##     - Actual LLM API calls (e.g., GPT-4, Claude) with per-query
##       outputs and uncertainty scores
##     - A real or quasi-real treatment-outcome dataset
##   None of this is available in this replication package.
##
## What this script demonstrates:
##   1. The full Algorithm 1 pipeline runs end-to-end
##   2. Cross-fitted pilot, water-filling allocation, integer
##      rounding, and stratified outcome model all work as
##      intended
##   3. The plug-in VOI estimator can underperform Unif at
##      n = 250 with realistic pilot noise (a finite-sample
##      limitation discussed in Section 8 of the manuscript)
## ============================================================

## Robust path resolution: works from repo root or scripts/
.calm_voi_path <- if (file.exists("R/calm_voi.R")) "R/calm_voi.R" else
                  if (file.exists("../R/calm_voi.R")) "../R/calm_voi.R" else
                  stop("Cannot find R/calm_voi.R. Run from repo root or scripts/.")
source(.calm_voi_path)


if (!dir.exists("outputs")) dir.create("outputs", recursive = TRUE)

set.seed(2026)

## ---- Synthetic trial DGP --------------------------------------
## 5 structured baseline features standing in for age, sex (centered),
## BMI, disease duration, baseline severity score
## 20-dimensional "latent text embedding" representing clinical notes
## A simple linear model maps embedding -> outcome with treatment effect
## and heterogeneous LLM noise driven by an "atypicality" score
## constructed from the embedding's L2 norm

simulate_trial <- function(n = 250) {
  # Structured features: standardized
  age   <- rnorm(n, 0, 1)
  sex   <- rbinom(n, 1, 0.55) - 0.55
  bmi   <- rnorm(n, 0, 1)
  durat <- rnorm(n, 0, 1)
  basev <- rnorm(n, 0, 1)
  X <- cbind(age, sex, bmi, durat, basev)
  colnames(X) <- c("age","sex","bmi","durat","basev")

  # Latent 20-d text embedding (highly informative about outcome)
  d_emb <- 20
  Z <- matrix(rnorm(n * d_emb), nrow = n)
  # Atypicality: high L2 norm => unusual notes => higher LLM uncertainty
  atypical <- (sqrt(rowSums(Z^2)) - sqrt(d_emb)) / sqrt(2 * d_emb)
  sigma_LLM_synth <- 0.6 + 2.5 * pnorm(atypical)  # 0.6 - 3.1 range

  # Treatment with 1:1 randomization
  A <- rbinom(n, 1, 0.5)
  pi_ <- 0.5

  # True outcome model (continuous symptom score; lower = better)
  beta_X    <- c(0.4, -0.3, 0.5, 0.2, 1.2)
  beta_Z    <- rnorm(d_emb, sd = 0.4)
  Xb        <- as.numeric(X %*% beta_X)
  Zsig      <- as.numeric(Z %*% beta_Z)
  Y0 <- 10 + Xb + Zsig + rnorm(n, 0, 1.0)
  Y1 <- 5  + Xb + Zsig + rnorm(n, 0, 1.0)  # treatment effect = -5
  Y  <- ifelse(A == 1, Y1, Y0)

  list(X = X, Z = Z, A = A, Y = Y, sigma_LLM_synth = sigma_LLM_synth)
}

## ---- Pseudo-LLM: pretrain a small MLP on Z -> outcome ---------
## In a real application this is replaced by an API call to a real LLM
pretrain_pseudo_llm <- function(n_train = 5000, d_emb = 20) {
  set.seed(99)
  Z_tr <- matrix(rnorm(n_train * d_emb), nrow = n_train)
  beta_Z_tr <- rnorm(d_emb, sd = 0.4)
  y_tr <- 10 + as.numeric(Z_tr %*% beta_Z_tr) + rnorm(n_train, 0, 0.5)
  ## Fit a smooth model: kernel ridge with random features
  d_rf <- 100
  W <- matrix(rnorm(d_emb * d_rf), nrow = d_emb)
  b <- runif(d_rf, 0, 2 * pi)
  phi_tr <- cos(Z_tr %*% W + matrix(b, nrow = n_train, ncol = d_rf,
                                    byrow = TRUE))
  ## Ridge regression
  lambda <- 1.0
  G <- t(phi_tr) %*% phi_tr + lambda * diag(d_rf)
  alpha <- solve(G, t(phi_tr) %*% y_tr)
  list(W = W, b = b, alpha = alpha,
       beta_Z_true = beta_Z_tr, d_rf = d_rf)
}

## Query pseudo-LLM on Z with temperature noise
## Returns predictions of shape (n, K_max)
query_pseudo_llm <- function(model, Z, K_max, sigma_LLM_synth) {
  n <- nrow(Z)
  phi <- cos(Z %*% model$W + matrix(model$b, nrow = n, ncol = model$d_rf,
                                     byrow = TRUE))
  ## Baseline (deterministic) prediction
  mu_LLM_baseline <- as.numeric(phi %*% model$alpha)
  ## Per-query noise: heteroscedastic Gaussian based on synthetic atypicality
  noise <- matrix(rnorm(n * K_max), nrow = n) * sigma_LLM_synth
  m_all <- mu_LLM_baseline + noise
  list(m_all = m_all, mu_baseline = mu_LLM_baseline)
}

## ---- Adapter to call our CALM-VOI machinery -------------------
run_one_trial <- function(trial, llm_model, K_budget) {
  X <- trial$X; A <- trial$A; Y <- trial$Y
  pi_ <- 0.5
  K_max_buffer <- 25
  qllm <- query_pseudo_llm(llm_model, trial$Z, K_max = K_max_buffer,
                           sigma_LLM_synth = trial$sigma_LLM_synth)
  m_all <- qllm$m_all
  m_first <- m_all[, 1]
  sigma_LLM <- trial$sigma_LLM_synth

  results <- list()

  ## AIPW
  r <- aipw_xfit(Y, A, X, M = NULL, pi_ = pi_, L = 5, include_m = FALSE)
  results$AIPW <- list(tau = r$tau, se = r$se, K = NA)

  ## CALM (K=1)
  r <- aipw_xfit(Y, A, X, M = m_first, pi_ = pi_, L = 5, include_m = TRUE)
  results$CALM <- list(tau = r$tau, se = r$se, K = rep(1L, length(Y)))

  ## Unif (K_i = K)
  K_u <- rep(as.integer(K_budget), length(Y))
  m_bar <- avg_predictions(m_all, K_u)
  r <- aipw_xfit(Y, A, X, M = m_bar, pi_ = pi_, L = 5, include_m = TRUE)
  results$Unif <- list(tau = r$tau, se = r$se, K = K_u)

  ## VOI plug-in
  pilot <- cross_fit_pilot_B(X, m_first, Y, A, pi_ = pi_, L = 5)
  q_hat <- pilot$B * sigma_LLM
  K_v_f <- water_filling(q_hat, length(Y) * K_budget, K_min = 1)
  K_v   <- integer_allocate(K_v_f, length(Y) * K_budget,
                            K_min = 1, K_max_cap = K_max_buffer)
  m_bar <- avg_predictions(m_all, K_v)
  r <- aipw_xfit(Y, A, X, M = m_bar, pi_ = pi_, L = 5, include_m = TRUE,
                 K_alloc = K_v, stratify_by_K = TRUE)
  results$VOI <- list(tau = r$tau, se = r$se, K = K_v)

  results
}

## ---- Main driver: single trial + mini-MC ----------------------
main_worked_example <- function() {
  cat("=== Pretraining pseudo-LLM... ===\n")
  llm_model <- pretrain_pseudo_llm()

  cat("=== Single illustrative trial ===\n")
  set.seed(1)
  trial <- simulate_trial(n = 250)
  res <- run_one_trial(trial, llm_model, K_budget = 5)
  tau_true <- -5
  cat(sprintf("True ATE: %.1f\n", tau_true))
  for (nm in names(res)) {
    tau <- res[[nm]]$tau; se <- res[[nm]]$se
    lo <- tau - 1.96 * se; hi <- tau + 1.96 * se
    Kmean <- if (is.na(res[[nm]]$K[1])) NA else mean(res[[nm]]$K)
    cat(sprintf("  %-6s tau=%.2f  SE=%.2f  CI=[%.2f, %.2f]  meanK=%s\n",
                nm, tau, se, lo, hi,
                if (is.na(Kmean)) "N/A" else sprintf("%.2f", Kmean)))
  }

  cat("\n=== Mini Monte Carlo (50 trials) ===\n")
  B_mini <- 50
  ests <- c("AIPW","CALM","Unif","VOI")
  tau_arr <- matrix(NA_real_, nrow = B_mini, ncol = length(ests))
  colnames(tau_arr) <- ests
  for (b in seq_len(B_mini)) {
    set.seed(100 + b)
    trial <- simulate_trial(n = 250)
    res <- run_one_trial(trial, llm_model, K_budget = 5)
    for (e in ests) tau_arr[b, e] <- res[[e]]$tau
  }
  cat(sprintf("%-8s %8s %8s %12s\n", "method", "mean", "SD", "vs AIPW (%)"))
  aipw_sd <- sd(tau_arr[, "AIPW"])
  for (e in ests) {
    m  <- mean(tau_arr[, e])
    s  <- sd(tau_arr[, e])
    red <- 100 * (1 - s^2 / aipw_sd^2)
    cat(sprintf("%-8s %8.2f %8.2f %12.1f\n", e, m, s, red))
  }

  ## Save the raw results for documentation
  write.csv(as.data.frame(tau_arr),
            "outputs/worked_example_mini_mc.csv", row.names = FALSE)
  cat("\nSaved worked_example_mini_mc.csv\n")
  cat("\nReminder: this is a SEMI-SYNTHETIC example.\n")
  cat("The 'LLM' is a random-features ridge model, not a real LLM.\n")
}

if (!interactive() && sys.nframe() == 0L) {
  main_worked_example()
}
