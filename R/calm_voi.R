## ============================================================
## calm_voi.R
## Core library for the CALM-VOI manuscript simulations.
## Implements:
##   - Three DGPs (reliability/leverage/cancellation)
##   - Parametric pilot and outcome model with x*m and m^2 terms
##   - Water-filling allocation with integer rounding
##   - Stratified cross-fitted AIPW estimator
##   - Five estimators (AIPW, CALM, Unif, Neyman, VOI plug-in) + VOI oracle
##
## Author: Claude (Anthropic), May 2026
## Notation matches the manuscript (CALM_VOI_manuscript.tex)
## ============================================================

suppressPackageStartupMessages({
  library(stats)
})

## ---- Globals ----------------------------------------------------
BETA_X <- c(1.0, 0.5, -0.3, 0.0, 0.0)
P_DIM  <- length(BETA_X)
TAU_TRUE <- 0.4

## ---- Baseline + treatment ---------------------------------------
generate_baseline <- function(n) {
  ## AR(1) covariance: Sigma_{jk} = 0.3^|j-k|
  Sigma <- 0.3 ^ abs(outer(1:P_DIM, 1:P_DIM, "-"))
  L <- chol(Sigma)
  X <- matrix(rnorm(n * P_DIM), nrow = n) %*% L
  Tstar <- rnorm(n)
  A <- rbinom(n, 1, 0.5)
  list(X = X, Tstar = Tstar, A = A)
}

## ---- Outcomes ---------------------------------------------------
## Returns Y, plus the TRUE per-patient leverages beta0_i, beta1_i.
## beta_ai = d mu*_a(X, mu^LLM) / d mu^LLM evaluated at the patient's
## mu^LLM_i. mu^LLM = 2*T* + X@beta_X, so T* = (mu^LLM - X@beta_X)/2.
##
## DGP 1: Y(a) = 0.4a + X@beta + 2*T* + eps; beta_0=beta_1=1
## DGP 2: Y(a) = 0.4a + X@beta + (2 + 0.5*X1)*T* + eps;
##              beta_0=beta_1=1 + 0.25*X1
## DGP 3: Y(a) = 0.4a + X@beta + [1 - 2a*I(X2<0)]*2*T* + eps
##              beta_0=1 always; beta_1=1 if X2>=0 else -1
make_outcomes <- function(dgp, X, Tstar, A) {
  n  <- nrow(X)
  eps <- rnorm(n, 0, sqrt(0.25))
  Xb <- as.numeric(X %*% BETA_X)
  if (identical(dgp, "1")) {
    Y0 <- Xb + 2 * Tstar + eps
    Y1 <- 0.4 + Xb + 2 * Tstar + eps
    beta0 <- rep(1, n); beta1 <- rep(1, n)
  } else if (identical(dgp, "2")) {
    coef <- 2 + 0.5 * X[, 1]
    Y0 <- Xb + coef * Tstar + eps
    Y1 <- 0.4 + Xb + coef * Tstar + eps
    beta0 <- 1 + 0.25 * X[, 1]; beta1 <- 1 + 0.25 * X[, 1]
  } else if (identical(dgp, "3")) {
    flip <- as.integer(X[, 2] < 0)
    coef0 <- rep(2, n)
    coef1 <- 2 * (1 - 2 * flip)              # -2 if X2<0
    Y0 <- Xb + coef0 * Tstar + eps
    Y1 <- 0.4 + Xb + coef1 * Tstar + eps
    beta0 <- rep(1, n)
    beta1 <- ifelse(flip == 1, -1, 1)
  } else stop("Unknown DGP: ", dgp)
  Y <- ifelse(A == 1, Y1, Y0)
  list(Y = Y, Y0 = Y0, Y1 = Y1, beta0 = beta0, beta1 = beta1)
}

## ---- LLM noise spec ---------------------------------------------
## spec: a number for homogeneous, "bimodal" for half 0.3 / half 2.0
make_sigma_eta <- function(spec, n) {
  if (identical(spec, "bimodal")) {
    h <- n %/% 2
    sig <- c(rep(0.3, h), rep(2.0, n - h))
    sig <- sample(sig)
    return(sig)
  }
  rep(as.numeric(spec), n)
}

## ---- LLM queries ------------------------------------------------
## m_i^(k) = 2*(T_i* + eta_i^(k)) + X_i@beta_X
##   so mu^LLM_i = 2*T_i* + X_i@beta_X and sigma_LLM_i = 2*sigma_eta_i
llm_queries <- function(X, Tstar, sigma_eta, K_max) {
  n <- nrow(X)
  Xb <- as.numeric(X %*% BETA_X)
  eta <- matrix(rnorm(n * K_max), nrow = n) * sigma_eta
  m <- 2 * (Tstar + eta) + Xb     # uses recycling via matrix; reset
  m <- 2 * (matrix(Tstar, nrow = n, ncol = K_max) + eta) +
       matrix(Xb, nrow = n, ncol = K_max)
  m
}

## ---- Parametric pilot / outcome model ---------------------------
## Form: mu(x, m) = theta0 + theta_x^T x + theta_m m
##                + theta_xm^T (x*m) + theta_mm m^2
design_matrix <- function(X, m) {
  ## [1, X (p), m, X*m (p), m^2]
  n <- nrow(X)
  Z <- cbind(1, X, m, X * m, m^2)
  Z
}

fit_arm_model <- function(X, m, Y) {
  Z <- design_matrix(X, m)
  ## solve(Z'Z) Z'Y but use qr.solve for stability
  qrZ <- qr(Z)
  qr.coef(qrZ, Y)
}

predict_value <- function(theta, X, m) {
  Z <- design_matrix(X, m)
  as.numeric(Z %*% theta)
}

predict_derivative <- function(theta, X, m) {
  ## d mu_hat / d m = theta_m + theta_xm^T X + 2*theta_mm*m
  p <- ncol(X)
  theta_m  <- theta[1 + p + 1]
  theta_xm <- theta[(1 + p + 2):(1 + 2 * p + 1)]
  theta_mm <- theta[1 + 2 * p + 2]
  as.numeric(theta_m + X %*% theta_xm + 2 * theta_mm * m)
}

## ---- Water-filling ----------------------------------------------
## Solves: K_i = max(K_min, c * q_i)  for q_i > 0
##         K_i = K_min                for q_i = 0
## subject to sum(K_i) = total_budget.
##
## Returns FLOAT allocations; use integer_allocate() to round.
water_filling <- function(q, total_budget, K_min = 1, tol = 1e-6,
                          max_iter = 200) {
  n <- length(q)
  K <- rep(K_min, n)
  pos <- q > 0
  n_pos <- sum(pos); n_zero <- n - n_pos
  if (n_pos == 0) return(rep(total_budget / n, n))
  q_pos <- q[pos]
  budget_pos <- total_budget - K_min * n_zero
  if (budget_pos <= K_min * n_pos) return(K)
  ## Binary search for c
  c_lo <- 0
  c_hi <- 1
  while (sum(pmax(K_min, c_hi * q_pos)) < budget_pos && c_hi < 1e12) {
    c_hi <- c_hi * 2
  }
  for (it in 1:max_iter) {
    c_mid <- 0.5 * (c_lo + c_hi)
    s <- sum(pmax(K_min, c_mid * q_pos))
    if (abs(s - budget_pos) < tol) break
    if (s < budget_pos) c_lo <- c_mid else c_hi <- c_mid
  }
  K[pos] <- pmax(K_min, c_mid * q_pos)
  K
}

integer_allocate <- function(K_float, total_budget, K_min = 1,
                             K_max_cap = NULL) {
  K <- pmax(K_min, floor(K_float))
  if (!is.null(K_max_cap)) K <- pmin(K, K_max_cap)
  diff_q <- total_budget - sum(K)
  if (diff_q > 0) {
    frac <- K_float - floor(K_float)
    eligible <- if (is.null(K_max_cap)) rep(TRUE, length(K)) else (K < K_max_cap)
    ord <- order(-frac)
    ord <- ord[eligible[ord]]
    if (length(ord) > 0) {
      take <- ord[seq_len(min(diff_q, length(ord)))]
      K[take] <- K[take] + 1
    }
  } else if (diff_q < 0) {
    excess <- -diff_q
    ord <- order(K_float - floor(K_float))
    removed <- 0
    for (i in ord) {
      if (removed >= excess) break
      if (K[i] > K_min) { K[i] <- K[i] - 1; removed <- removed + 1 }
    }
    while (sum(K) > total_budget) {
      i <- which.max(K)
      if (K[i] > K_min) K[i] <- K[i] - 1 else break
    }
  }
  K
}

## ---- Cross-fitting helpers --------------------------------------
make_folds <- function(n, L = 5) {
  ## Random fold assignment with sizes as equal as possible
  fold_id <- ((sample.int(n) - 1L) %% L) + 1L
  fold_id
}

## Compute averaged covariate based on per-patient allocation.
## Vectorized: builds a 0/1 mask matrix and divides row-sums by K_alloc.
avg_predictions <- function(m_all, K_alloc) {
  n <- nrow(m_all); K_max <- ncol(m_all)
  ## mask[i, k] = 1 iff k <= K_alloc[i]
  mask <- outer(K_alloc, seq_len(K_max), FUN = ">=")
  storage.mode(mask) <- "double"
  rowSums(m_all * mask) / K_alloc
}

## ---- Cross-fitted AIPW with optional stratified outcome model ----
## stratify_by_K: if TRUE, fit separate outcome models for
##                strata defined by (K_alloc==1) vs (K_alloc>1).
## include_m: if FALSE, fit linear-in-X outcome model (used by AIPW).
aipw_xfit <- function(Y, A, X, M, pi_ = 0.5, L = 5,
                      include_m = TRUE, K_alloc = NULL,
                      stratify_by_K = FALSE) {
  n <- length(Y)
  fold <- make_folds(n, L)
  if (stratify_by_K && !is.null(K_alloc)) {
    stratum <- as.integer(K_alloc > 1)  # 0 = "K_i=1", 1 = "K_i>1"
  } else {
    stratum <- rep(0L, n)
  }
  mu1 <- numeric(n); mu0 <- numeric(n)
  min_per_stratum <- 2 * P_DIM + 4
  for (ell in seq_len(L)) {
    te <- which(fold == ell)
    tr <- which(fold != ell)
    for (a in c(0L, 1L)) {
      mu_out <- if (a == 1L) "mu1" else "mu0"
      for (s in unique(stratum)) {
        tr_mask <- tr[A[tr] == a & stratum[tr] == s]
        te_mask <- te[stratum[te] == s]
        if (length(tr_mask) < min_per_stratum) {
          tr_mask <- tr[A[tr] == a]
        }
        if (length(te_mask) == 0L) next
        if (include_m) {
          theta <- fit_arm_model(X[tr_mask, , drop = FALSE],
                                 M[tr_mask], Y[tr_mask])
          pred  <- predict_value(theta, X[te_mask, , drop = FALSE],
                                 M[te_mask])
        } else {
          ## linear-in-X for AIPW baseline
          Z_tr <- cbind(1, X[tr_mask, , drop = FALSE])
          theta_lin <- qr.coef(qr(Z_tr), Y[tr_mask])
          Z_te <- cbind(1, X[te_mask, , drop = FALSE])
          pred <- as.numeric(Z_te %*% theta_lin)
        }
        if (a == 1L) mu1[te_mask] <- pred else mu0[te_mask] <- pred
      }
    }
  }
  scores <- A * (Y - mu1) / pi_ -
            (1 - A) * (Y - mu0) / (1 - pi_) +
            mu1 - mu0
  tau <- mean(scores)
  se  <- sqrt(var(scores) / n)
  list(tau = tau, se = se, scores = scores)
}

## ---- Pilot for plug-in B_i --------------------------------------
cross_fit_pilot_B <- function(X, m_first, Y, A, pi_ = 0.5, L = 5) {
  n <- length(Y)
  fold <- make_folds(n, L)
  beta0_hat <- numeric(n); beta1_hat <- numeric(n)
  for (ell in seq_len(L)) {
    te <- which(fold == ell)
    tr <- which(fold != ell)
    for (a in c(0L, 1L)) {
      mask <- tr[A[tr] == a]
      theta <- fit_arm_model(X[mask, , drop = FALSE],
                             m_first[mask], Y[mask])
      pred  <- predict_derivative(theta, X[te, , drop = FALSE],
                                  m_first[te])
      if (a == 1L) beta1_hat[te] <- pred else beta0_hat[te] <- pred
    }
  }
  B_hat <- abs((1 - pi_) * beta1_hat + pi_ * beta0_hat) /
           sqrt(pi_ * (1 - pi_))
  list(B = B_hat, beta0 = beta0_hat, beta1 = beta1_hat)
}

## ---- One Monte Carlo replication --------------------------------
## Returns a named list of estimator results.
## sigma_for_alloc: which sigma to use for allocation. Options:
##   "true"  : use true sigma_LLM (oracle uncertainty, default)
##   numeric vector of length n: pre-computed plug-in or noised sigma
## winsor_q: optional quantile for winsorizing q_hat
one_replication <- function(n, K_budget, dgp, sigma_eta_spec,
                            L = 5, K_max_buffer = 25,
                            sigma_for_alloc = "true",
                            winsor_q = NULL) {
  base   <- generate_baseline(n)
  X      <- base$X; Tstar <- base$Tstar; A <- base$A
  sigma_eta <- make_sigma_eta(sigma_eta_spec, n)
  out    <- make_outcomes(dgp, X, Tstar, A)
  Y <- out$Y; beta0_true <- out$beta0; beta1_true <- out$beta1
  pi_ <- 0.5
  sigma_LLM_true <- 2 * sigma_eta
  B_true <- abs((1 - pi_) * beta1_true + pi_ * beta0_true) /
            sqrt(pi_ * (1 - pi_))
  q_true <- B_true * sigma_LLM_true

  m_all <- llm_queries(X, Tstar, sigma_eta, K_max_buffer)
  m_first <- m_all[, 1]

  if (identical(sigma_for_alloc, "true")) {
    sigma_alloc <- sigma_LLM_true
  } else {
    sigma_alloc <- sigma_for_alloc
  }

  results <- list()

  ## 1. AIPW
  r <- aipw_xfit(Y, A, X, M = NULL, pi_ = pi_, L = L, include_m = FALSE)
  results$AIPW <- list(tau = r$tau, se = r$se, K = NA)

  ## 2. CALM (K=1)
  r <- aipw_xfit(Y, A, X, M = m_first, pi_ = pi_, L = L,
                 include_m = TRUE)
  results$CALM <- list(tau = r$tau, se = r$se, K = rep(1L, n))

  ## 3. Unif (K_i = K)
  K_unif <- rep(as.integer(K_budget), n)
  m_bar  <- avg_predictions(m_all, K_unif)
  r <- aipw_xfit(Y, A, X, M = m_bar, pi_ = pi_, L = L, include_m = TRUE)
  results$Unif <- list(tau = r$tau, se = r$se, K = K_unif)

  ## 4. Neyman (K_i propto sigma_alloc)
  K_n_f  <- water_filling(sigma_alloc, n * K_budget, K_min = 1)
  K_n    <- integer_allocate(K_n_f, n * K_budget, K_min = 1,
                             K_max_cap = K_max_buffer)
  m_bar  <- avg_predictions(m_all, K_n)
  r <- aipw_xfit(Y, A, X, M = m_bar, pi_ = pi_, L = L, include_m = TRUE,
                 K_alloc = K_n, stratify_by_K = TRUE)
  results$Neyman <- list(tau = r$tau, se = r$se, K = K_n)

  ## 5. VOI plug-in (K_i propto B_hat * sigma_alloc)
  pilot   <- cross_fit_pilot_B(X, m_first, Y, A, pi_ = pi_, L = L)
  q_hat   <- pilot$B * sigma_alloc
  if (!is.null(winsor_q)) {
    cap <- quantile(q_hat, winsor_q, na.rm = TRUE)
    q_hat <- pmin(q_hat, cap)
  }
  K_v_f   <- water_filling(q_hat, n * K_budget, K_min = 1)
  K_v     <- integer_allocate(K_v_f, n * K_budget, K_min = 1,
                              K_max_cap = K_max_buffer)
  m_bar   <- avg_predictions(m_all, K_v)
  r <- aipw_xfit(Y, A, X, M = m_bar, pi_ = pi_, L = L, include_m = TRUE,
                 K_alloc = K_v, stratify_by_K = TRUE)
  results$VOI <- list(tau = r$tau, se = r$se, K = K_v, B_hat = pilot$B)

  ## 5b. VOI oracle (true B_i and true sigma_LLM)
  K_v_or_f <- water_filling(q_true, n * K_budget, K_min = 1)
  K_v_or   <- integer_allocate(K_v_or_f, n * K_budget, K_min = 1,
                               K_max_cap = K_max_buffer)
  m_bar    <- avg_predictions(m_all, K_v_or)
  r <- aipw_xfit(Y, A, X, M = m_bar, pi_ = pi_, L = L, include_m = TRUE,
                 K_alloc = K_v_or, stratify_by_K = TRUE)
  results$VOI_oracle <- list(tau = r$tau, se = r$se, K = K_v_or,
                             B_true = B_true)

  attr(results, "X") <- X
  attr(results, "B_true") <- B_true
  attr(results, "sigma_LLM_true") <- sigma_LLM_true
  attr(results, "B_hat") <- pilot$B
  results
}
