## ============================================================
## calm_voi.R
## Core library for the CALM-VOI manuscript simulations.
## Implements:
##   - Three DGPs (reliability/leverage/cancellation)
##   - Pooled parametric pilot and outcome model with x*m and m^2 terms
##   - Exact greedy integer allocation within analysis folds
##   - Shrinkage allocation and the cross-fitted stability diagnostic
##   - Stratified cross-fitted AIPW estimator
##   - AIPW, CALM, Unif, Neyman, VOI plug-in, VOI oracle, VOI
##     winsorized, VOI diagnostic-gated, VOI shrink
##
## Notation matches the manuscript.
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
## spec: a number for homogeneous, "bimodal" for half 0.3 / half 2.0,
##       "lognormal" for continuous heterogeneity.
make_sigma_eta <- function(spec, n) {
  if (identical(spec, "bimodal")) {
    h <- n %/% 2
    sig <- c(rep(0.3, h), rep(2.0, n - h))
    sig <- sample(sig)
    return(sig)
  }
  if (identical(spec, "lognormal")) {
    ## Continuous heterogeneity: median 1.0, CV ~ 0.53, so comparable
    ## in spread to "bimodal" but with genuine rank structure.
    ## Required by the rank-corruption regimes in run_sensitivity.R,
    ## which are a no-op when sigma is constant and near-degenerate
    ## when sigma takes only two distinct values.
    return(exp(rnorm(n, 0, 0.5)))
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

derivative_matrix <- function(X, m) {
  ## Derivative of design_matrix(X, m) with respect to m.
  n <- nrow(X)
  cbind(0, matrix(0, nrow = n, ncol = ncol(X)), 1, X, 2 * m)
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

## Joint pilot for arm-specific surrogate derivatives. Centering treatment at
## the known randomization probability separates the shared prognostic curve
## from arm differences. The shared curve uses every training observation;
## ridge regularization partially pools the arm differences. The first column
## of Z must be an intercept. interaction_cols index the nonintercept columns
## of Z that may differ by treatment.
fit_pooled_arm_model <- function(Z, Y, A, pi_ = 0.5,
                                 interaction_cols = NULL,
                                 ridge_main = 1e-6,
                                 ridge_interaction = NULL) {
  Z <- as.matrix(Z)
  n <- nrow(Z)
  p <- ncol(Z) - 1L
  if (n < 2L || p < 1L || length(Y) != n || length(A) != n ||
      any(!is.finite(Z)) || any(!is.finite(Y)) ||
      any(!A %in% c(0, 1)) || !is.finite(pi_) || pi_ <= 0 || pi_ >= 1 ||
      any(abs(Z[, 1L] - 1) > 1e-10)) {
    stop("Invalid inputs to fit_pooled_arm_model")
  }
  if (is.null(interaction_cols)) interaction_cols <- seq_len(p)
  interaction_cols <- sort(unique(as.integer(interaction_cols)))
  if (anyNA(interaction_cols) || !length(interaction_cols) ||
      any(interaction_cols < 1L | interaction_cols > p)) {
    stop("interaction_cols must index nonintercept columns of Z")
  }
  if (!is.finite(ridge_main) || ridge_main < 0) {
    stop("ridge_main must be finite and nonnegative")
  }
  if (is.null(ridge_interaction)) {
    ridge_interaction <- 0.5 / sqrt(n)
  }
  if (!is.finite(ridge_interaction) || ridge_interaction < 0) {
    stop("ridge_interaction must be finite and nonnegative")
  }

  U <- Z[, -1L, drop = FALSE]
  center <- colMeans(U)
  U_centered <- sweep(U, 2L, center, "-")
  scale <- sqrt(colMeans(U_centered^2))
  scale[!is.finite(scale) | scale < sqrt(.Machine$double.eps)] <- 1
  U_scaled <- sweep(U_centered, 2L, scale, "/")
  Ac <- A - pi_
  Ui <- U_scaled[, interaction_cols, drop = FALSE]
  W <- cbind(1, U_scaled, Ac, Ac * Ui)

  ## Penalize standardized slopes, but not the marginal or treatment
  ## intercept. The n multiplier expresses ridge on the mean squared-error
  ## scale and makes the defaults comparable across sample sizes.
  penalty <- c(
    0, rep(ridge_main, p), 0,
    rep(ridge_interaction, length(interaction_cols))
  )
  augmented_W <- rbind(W, diag(sqrt(n * penalty), nrow = length(penalty)))
  augmented_Y <- c(Y, rep(0, length(penalty)))
  coef <- qr.coef(qr(augmented_W), augmented_Y)
  if (anyNA(coef)) {
    stop("Rank-deficient pooled pilot fit; remove redundant unpenalized terms")
  }
  list(
    shared = coef[seq_len(p + 1L)],
    difference = coef[(p + 2L):length(coef)],
    center = center,
    scale = scale,
    interaction_cols = interaction_cols,
    pi = pi_,
    ridge_main = ridge_main,
    ridge_interaction = ridge_interaction
  )
}

predict_pooled_derivatives <- function(fit, dZ) {
  dZ <- as.matrix(dZ)
  p <- length(fit$scale)
  if (ncol(dZ) != p + 1L || any(!is.finite(dZ))) {
    stop("dZ is incompatible with the pooled pilot fit")
  }
  dU <- sweep(dZ[, -1L, drop = FALSE], 2L, fit$scale, "/")
  shared <- as.numeric(cbind(dZ[, 1L], dU) %*% fit$shared)
  dUi <- dU[, fit$interaction_cols, drop = FALSE]
  difference <- as.numeric(
    cbind(dZ[, 1L], dUi) %*% fit$difference
  )
  cbind(
    beta0 = shared - fit$pi * difference,
    beta1 = shared + (1 - fit$pi) * difference
  )
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
  if (!n || any(!is.finite(q)) || any(q < 0)) {
    stop("q must be a nonempty finite nonnegative vector")
  }
  if (length(total_budget) != 1L || !is.finite(total_budget) ||
      total_budget < n * K_min) {
    stop("total_budget must be finite and at least length(q) * K_min")
  }
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
  n <- length(K_float)
  if (!n || any(!is.finite(K_float))) {
    stop("K_float must be a nonempty finite vector")
  }
  if (length(total_budget) != 1L || !is.finite(total_budget) ||
      total_budget != round(total_budget)) {
    stop("total_budget must be a finite integer")
  }
  if (length(K_min) != 1L || !is.finite(K_min) ||
      K_min < 1 || K_min != round(K_min)) {
    stop("K_min must be a positive integer")
  }
  if (!is.null(K_max_cap) &&
      (length(K_max_cap) != 1L || K_max_cap < K_min ||
       (is.finite(K_max_cap) && K_max_cap != round(K_max_cap)))) {
    stop("K_max_cap must be NULL, Inf, or an integer no smaller than K_min")
  }
  if (total_budget < n * K_min ||
      (!is.null(K_max_cap) && is.finite(K_max_cap) &&
       total_budget > n * K_max_cap)) {
    stop("total_budget is infeasible under K_min and K_max_cap")
  }
  K <- pmax(K_min, floor(K_float))
  if (!is.null(K_max_cap)) K <- pmin(K, K_max_cap)
  diff_q <- total_budget - sum(K)
  if (diff_q > 0) {
    frac <- K_float - floor(K_float)
    ord_all <- order(-frac)
    ## Repeated passes. A single pass adds at most one query per
    ## eligible patient, which silently leaves the budget unspent
    ## whenever K_max_cap binds for many patients (this happens under
    ## extreme allocation scores: at sigma-noise sd 1.5 the shortfall
    ## averaged 6% of nK and reached 35% in the worst replication).
    ## Terminates when the budget is met or every patient is capped.
    repeat {
      if (diff_q <= 0) break
      eligible <- if (is.null(K_max_cap)) rep(TRUE, length(K)) else (K < K_max_cap)
      ord <- ord_all[eligible[ord_all]]
      if (length(ord) == 0L) break
      take <- ord[seq_len(min(diff_q, length(ord)))]
      K[take] <- K[take] + 1
      diff_q <- total_budget - sum(K)
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
  if (sum(K) != total_budget) {
    stop("integer allocation could not satisfy the requested budget and caps")
  }
  as.integer(K)
}

## ---- Exact integer allocation by greedy marginal gain ------------
## Solves exactly
##   min sum_i q_i^2 / K_i   s.t. sum_i K_i = total_budget,
##                                K_min <= K_i <= K_max, K_i integer.
## Starting from K_i = K_min, each further query goes to the patient
## with the largest marginal reduction
##   Delta_i(k) = q_i^2/k - q_i^2/(k+1) = q_i^2 / (k(k+1)).
## Delta_i(k) is decreasing in k, so the greedy order is optimal for
## this separable convex program: no rounding step is needed.
##
## Rather than performing total_budget heap pops, we threshold. Patient
## i takes every step whose gain is at least t, i.e. all k with
## k(k+1) <= q_i^2/t, which has the closed form
##   k <= ke_i(t) = (-1 + sqrt(1 + 4 q_i^2 / t)) / 2,
## so K_i(t) = K_min + max(0, floor(ke_i(t)) - K_min + 1). K_i(t) is
## non-increasing in t; a geometric bisection lands at the budget up to
## boundary ties, which are then settled by explicit marginal steps.
greedy_allocate <- function(q, total_budget, K_min = 1L, K_max = Inf) {
  n  <- length(q)
  if (!n || any(!is.finite(q)) || any(q < 0)) {
    stop("q must be a nonempty finite nonnegative vector")
  }
  if (length(K_min) != 1L || !is.finite(K_min) ||
      K_min < 1 || K_min != as.integer(K_min)) {
    stop("K_min must be a positive integer")
  }
  if (length(K_max) != 1L || K_max < K_min ||
      (is.finite(K_max) && K_max != as.integer(K_max))) {
    stop("K_max must be Inf or an integer no smaller than K_min")
  }
  if (length(total_budget) != 1L || !is.finite(total_budget) ||
      total_budget != round(total_budget)) {
    stop("total_budget must be a finite integer")
  }
  ## The optimizer is invariant to a common positive rescaling. Scale
  ## before squaring to avoid overflow for otherwise finite scores.
  q_num <- as.numeric(q)
  q_scale <- max(q_num)
  q2 <- if (q_scale > 0) (q_num / q_scale)^2 else numeric(n)
  total_budget <- as.integer(total_budget)
  if (total_budget < n * K_min ||
      (is.finite(K_max) && total_budget > n * K_max)) {
    stop("total_budget is infeasible under K_min and K_max")
  }
  if (total_budget == n * K_min) return(rep(as.integer(K_min), n))
  if (all(q2 <= 0)) {
    ## Objective is zero for any allocation, but the budget is still
    ## spent in a real trial; spread it as evenly as the cap allows
    ## rather than silently returning n*K_min queries.
    base <- min(floor(total_budget / n), K_max)
    K    <- rep(base, n)
    rem  <- total_budget - sum(K)
    if (rem > 0) {
      add <- which(K < K_max)[seq_len(min(rem, sum(K < K_max)))]
      K[add] <- K[add] + 1
    }
    return(as.integer(K))
  }
  K_at <- function(t) {
    ke <- numeric(n)
    pos <- q2 > 0
    ke[pos] <- floor((-1 + sqrt(1 + 4 * q2[pos] / t)) / 2)
    pmin(K_min + pmax(0, ke - K_min + 1), K_max)
  }
  ## bracket the threshold
  t_hi <- max(q2)
  while (sum(K_at(t_hi)) > total_budget && is.finite(t_hi)) t_hi <- t_hi * 2
  t_lo <- t_hi
  while (sum(K_at(t_lo)) < total_budget && t_lo > 1e-300) t_lo <- t_lo / 2
  for (it in seq_len(200)) {
    tm <- sqrt(t_lo) * sqrt(t_hi)
    if (!is.finite(tm) || tm <= 0) break
    if (sum(K_at(tm)) >= total_budget) t_lo <- tm else t_hi <- tm
  }
  K <- K_at(t_lo)                     # sum(K) >= total_budget
  ## settle ties: drop the least valuable assigned steps
  guard <- 0L
  while (sum(K) > total_budget && guard < 10L * n + 200L) {
    guard <- guard + 1L
    elig <- which(K > K_min)
    if (!length(elig)) break
    j <- elig[which.min(q2[elig] / ((K[elig] - 1) * K[elig]))]
    K[j] <- K[j] - 1
  }
  ## if the cap left the budget short, add the most valuable steps
  guard <- 0L
  while (sum(K) < total_budget && guard < 10L * n + 200L) {
    guard <- guard + 1L
    elig <- which(K < K_max)
    if (!length(elig)) break
    j <- elig[which.max(q2[elig] / (K[elig] * (K[elig] + 1)))]
    K[j] <- K[j] + 1
  }
  if (sum(K) != total_budget) {
    stop("greedy allocation failed to spend the exact feasible budget")
  }
  as.integer(K)
}

## Allocate separately within prespecified analysis folds.  This prevents
## the budget threshold for a held-out patient from depending on allocation
## scores in other folds whose pilot fits may have used that patient's
## outcome.  Each fold receives exactly length(fold) * K_budget queries, so
## the grand total remains n * K_budget.
greedy_allocate_by_fold <- function(q, fold_id, K_budget,
                                    K_min = 1L, K_max = Inf) {
  stopifnot(length(q) == length(fold_id))
  K <- integer(length(q))
  for (ell in sort(unique(fold_id))) {
    ii <- which(fold_id == ell)
    K[ii] <- greedy_allocate(q[ii], length(ii) * K_budget,
                             K_min = K_min, K_max = K_max)
  }
  K
}

## ---- Shrinkage allocation (Proposition: shrinkage regret) --------
## Continuous shrinkage target interpolating uniform and plug-in:
##   Ktilde_i(lambda) = lambda * Khat_i^c + (1 - lambda) * K,
## with Khat_i^c = max(1, chat * qhat_i) the continuous water-filling
## solution. Because sum_i Khat_i^c = nK, the target also sums to nK
## for every lambda, so it is budget-feasible throughout.
##
## We interpolate between the exact integer plug-in allocation and uniform,
## then use largest-remainder integerization.  This preserves the endpoints
## exactly: lambda=0 is uniform and lambda=1 is the exact plug-in allocation.
shrinkage_allocate <- function(q_hat, n, K_budget, lambda,
                               K_min = 1L, K_max = Inf) {
  total <- n * K_budget
  K_plugin <- greedy_allocate(q_hat, total, K_min = K_min, K_max = K_max)
  target <- lambda * K_plugin + (1 - lambda) * K_budget
  integer_allocate(target, total, K_min = K_min, K_max_cap = K_max)
}

shrinkage_allocate_by_fold <- function(q_hat, fold_id, K_budget, lambda,
                                       K_min = 1L, K_max = Inf) {
  stopifnot(length(q_hat) == length(fold_id))
  K <- integer(length(q_hat))
  for (ell in sort(unique(fold_id))) {
    ii <- which(fold_id == ell)
    K[ii] <- shrinkage_allocate(q_hat[ii], length(ii), K_budget, lambda,
                                K_min = K_min, K_max = K_max)
  }
  K
}

## Apply a score cap separately within analysis folds. A single global
## empirical quantile would let a training-fold score change the cap in
## another fold and reintroduce the reverse-leakage path that fold-local
## budgeting is designed to remove.
winsorize_by_fold <- function(q, fold_id, prob = 0.95) {
  stopifnot(length(q) == length(fold_id), prob > 0, prob <= 1)
  out <- as.numeric(q)
  for (ell in sort(unique(fold_id))) {
    ii <- which(fold_id == ell)
    cap <- quantile(out[ii], prob, na.rm = TRUE, names = FALSE)
    out[ii] <- pmin(out[ii], cap)
  }
  out
}

## ---- Shrinkage weight from the stability diagnostic --------------
## The manuscript specifies that lambda is set by the pre-specified
## cross-fitted stability diagnostic but does not give the map from the
## diagnostic to lambda in [0,1]; the rule below is our construction.
## It is linear in the dispersion statistic and anchored at the
## supplement's pre-specified fallback threshold: lambda = 1 (full
## plug-in) when CV90 <= cv_lo, lambda = 0 (full fallback) when
## CV90 >= cv_hi = 0.5, interpolating in between. Within the shrinkage
## path, lambda = 0 is exactly uniform; the simulation wrapper may
## separately route that case to Neyman when uncertainty is judged
## heterogeneous and reliable.
select_lambda <- function(CV90, cv_lo = 0.2, cv_hi = 0.5) {
  if (!is.finite(CV90)) return(0)
  max(0, min(1, (cv_hi - CV90) / (cv_hi - cv_lo)))
}

## ---- Cross-fitting helpers --------------------------------------
make_folds <- function(n, L = 5) {
  ## Random fold assignment with sizes as equal as possible
  fold_id <- ((sample.int(n) - 1L) %% L) + 1L
  fold_id
}

## Empirical fourth-moment Monte Carlo standard error for a sample
## variance.  Unlike the normal-theory shortcut s^2*sqrt(2/(B-1)),
## this remains informative when the Monte Carlo distribution is
## heavy-tailed.
mcse_sample_variance <- function(x) {
  x <- x[is.finite(x)]
  B <- length(x)
  if (B < 4L) return(NA_real_)
  s2 <- var(x)
  m4 <- mean((x - mean(x))^4)
  sqrt(max(0, (m4 - ((B - 3) / (B - 1)) * s2^2) / B))
}

## Compute averaged covariate based on per-patient allocation.
## Vectorized: builds a 0/1 mask matrix and divides row-sums by K_alloc.
avg_predictions <- function(m_all, K_alloc) {
  n <- nrow(m_all); K_max <- ncol(m_all)
  if (length(K_alloc) != n || any(!is.finite(K_alloc)) ||
      any(K_alloc < 1) || any(K_alloc > K_max) ||
      any(K_alloc != round(K_alloc))) {
    stop("K_alloc must contain one feasible positive integer per row of m_all")
  }
  ## mask[i, k] = 1 iff k <= K_alloc[i]
  mask <- outer(K_alloc, seq_len(K_max), FUN = ">=")
  storage.mode(mask) <- "double"
  rowSums(m_all * mask) / K_alloc
}

## Primary frozen-score analysis from Algorithm 1. The caller must
## construct q_frozen on independent design data before generating the
## analysis query matrix m_all. This function applies one cohort-wide
## exact integer budget and never estimates an allocation score from Y,
## A, or m_all.
frozen_voi_analysis <- function(Y, A, X, m_all, q_frozen, K_budget,
                                K_max = ncol(m_all), pi_ = 0.5, L = 5,
                                fold_id = NULL, stratify_by_K = TRUE) {
  n <- length(Y)
  if (length(A) != n || nrow(X) != n || nrow(m_all) != n ||
      length(q_frozen) != n || any(!is.finite(q_frozen)) ||
      any(q_frozen < 0)) {
    stop("Frozen analysis inputs have incompatible dimensions or scores")
  }
  K_max <- min(as.integer(K_max), ncol(m_all))
  total_budget <- n * K_budget
  if (!is.finite(total_budget) || total_budget != round(total_budget) ||
      total_budget < n || total_budget > n * K_max) {
    stop("n * K_budget must be a feasible integer total")
  }
  K_alloc <- greedy_allocate(q_frozen, as.integer(total_budget),
                             K_min = 1L, K_max = K_max)
  M <- avg_predictions(m_all, K_alloc)
  fit <- aipw_xfit(Y, A, X, M, pi_ = pi_, L = L, include_m = TRUE,
                   K_alloc = K_alloc, stratify_by_K = stratify_by_K,
                   fold_id = fold_id)
  c(fit, list(K_alloc = K_alloc, M = M, q_frozen = q_frozen))
}

## ---- Cross-fitted AIPW with optional stratified outcome model ----
## stratify_by_K: if TRUE, fit separate outcome models for
##                strata defined by (K_alloc==1) vs (K_alloc>1).
## include_m: if FALSE, fit linear-in-X outcome model (used by AIPW).
aipw_xfit <- function(Y, A, X, M, pi_ = 0.5, L = 5,
                      include_m = TRUE, K_alloc = NULL,
                      stratify_by_K = FALSE, fold_id = NULL) {
  n <- length(Y)
  fold <- if (is.null(fold_id)) make_folds(n, L) else as.integer(fold_id)
  if (length(fold) != n || anyNA(fold)) stop("Invalid fold_id")
  L <- length(unique(fold))
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
          ## qr.coef returns NA for aliased columns, i.e. exactly when
          ## the training block is rank deficient. The size guard above
          ## is only one observation above the 2p+3 column count, so a
          ## thin stratum can satisfy it and still be deficient (seen
          ## in worked_example.R: a 14-row block of rank 12/13). Left
          ## unhandled the NA propagates through mu into tau, and a
          ## single NA replication makes the whole Monte Carlo column
          ## NA. Retry pooled across strata within the arm, then fall
          ## back to a linear-in-X fit.
          if (anyNA(theta)) {
            tr_pool <- tr[A[tr] == a]
            theta <- fit_arm_model(X[tr_pool, , drop = FALSE],
                                   M[tr_pool], Y[tr_pool])
          }
          if (anyNA(theta)) {
            ## The nonlinear pooled fit has still failed, so use the
            ## pooled within-arm sample for the simpler linear fallback
            ## as well. Reusing the original thin stratum here can leave
            ## even the linear design underdetermined.
            Z_tr  <- cbind(1, X[tr_pool, , drop = FALSE])
            th_l  <- qr.coef(qr(Z_tr), Y[tr_pool])
            th_l[is.na(th_l)] <- 0
            pred  <- as.numeric(cbind(1, X[te_mask, , drop = FALSE]) %*% th_l)
          } else {
            pred <- predict_value(theta, X[te_mask, , drop = FALSE],
                                  M[te_mask])
          }
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
cross_fit_pilot_B <- function(X, m_first, Y, A, pi_ = 0.5, L = 5,
                              fold_id = NULL,
                              pilot = c("pooled", "separate"),
                              ridge_main = 1e-6,
                              ridge_interaction = NULL) {
  n <- length(Y)
  pilot <- match.arg(pilot)
  fold <- if (is.null(fold_id)) make_folds(n, L) else as.integer(fold_id)
  if (length(fold) != n || anyNA(fold)) stop("Invalid fold_id")
  L <- length(unique(fold))
  beta0_hat <- numeric(n); beta1_hat <- numeric(n)
  for (ell in seq_len(L)) {
    te <- which(fold == ell)
    tr <- which(fold != ell)
    if (pilot == "pooled") {
      fit <- fit_pooled_arm_model(
        design_matrix(X[tr, , drop = FALSE], m_first[tr]),
        Y[tr], A[tr], pi_ = pi_, ridge_main = ridge_main,
        ridge_interaction = ridge_interaction
      )
      beta <- predict_pooled_derivatives(
        fit, derivative_matrix(X[te, , drop = FALSE], m_first[te])
      )
      beta0_hat[te] <- beta[, "beta0"]
      beta1_hat[te] <- beta[, "beta1"]
    } else {
      for (a in c(0L, 1L)) {
        mask <- tr[A[tr] == a]
        theta <- fit_arm_model(X[mask, , drop = FALSE],
                               m_first[mask], Y[mask])
        if (anyNA(theta)) {
          stop("Rank-deficient pilot fit; increase the pilot sample or simplify ",
               "the pilot basis")
        }
        pred <- predict_derivative(theta, X[te, , drop = FALSE],
                                   m_first[te])
        if (a == 1L) beta1_hat[te] <- pred else beta0_hat[te] <- pred
      }
    }
  }
  B_hat <- abs((1 - pi_) * beta1_hat + pi_ * beta0_hat) /
           sqrt(pi_ * (1 - pi_))
  list(B = B_hat, beta0 = beta0_hat, beta1 = beta1_hat)
}

## ---- Pre-specified pilot-quality diagnostic ----------------------
## Implements Supplement "Pre-Specified Diagnostic for Plug-in Pilot
## Quality": refit the pilot with L_diag = 10 delete-fold fits and
## evaluate every fit at every patient, giving ten leverage estimates
## B_hat_i^(l). Exactly one is out-of-fold for a given patient; nine
## training sets include that patient. The patient-level dispersion is
##   CV_i  = sd(B_hat_i^(l)) / |mean(B_hat_i^(l))|
## and the reported statistic is CV90 = quantile_0.90(CV_i).
pilot_B_dispersion <- function(X, m_first, Y, A, pi_ = 0.5, L_diag = 10,
                               pilot = c("pooled", "separate"),
                               ridge_main = 1e-6,
                               ridge_interaction = NULL) {
  n <- length(Y)
  pilot <- match.arg(pilot)
  fold <- make_folds(n, L_diag)
  Bmat <- matrix(NA_real_, nrow = n, ncol = L_diag)
  min_fit <- 2 * P_DIM + 4
  for (ell in seq_len(L_diag)) {
    tr <- which(fold != ell)
    if (pilot == "pooled") {
      fit <- fit_pooled_arm_model(
        design_matrix(X[tr, , drop = FALSE], m_first[tr]),
        Y[tr], A[tr], pi_ = pi_, ridge_main = ridge_main,
        ridge_interaction = ridge_interaction
      )
      beta <- predict_pooled_derivatives(
        fit, derivative_matrix(X, m_first)
      )
      Bmat[, ell] <- abs(
        (1 - pi_) * beta[, "beta1"] + pi_ * beta[, "beta0"]
      ) / sqrt(pi_ * (1 - pi_))
    } else {
      b <- list()
      ok <- TRUE
      for (a in c(0L, 1L)) {
        mask <- tr[A[tr] == a]
        if (length(mask) < min_fit) { ok <- FALSE; break }
        theta <- fit_arm_model(
          X[mask, , drop = FALSE], m_first[mask], Y[mask]
        )
        if (anyNA(theta)) { ok <- FALSE; break }
        ## evaluate this fold's coefficients at all patients
        b[[as.character(a)]] <- predict_derivative(theta, X, m_first)
      }
      if (!ok) next
      Bmat[, ell] <- abs((1 - pi_) * b[["1"]] + pi_ * b[["0"]]) /
                     sqrt(pi_ * (1 - pi_))
    }
  }
  mu  <- rowMeans(Bmat, na.rm = TRUE)
  sdv <- apply(Bmat, 1, sd, na.rm = TRUE)
  ## CV is undefined where the mean leverage is ~0. Those patients are
  ## exactly the arm-cancellation cases (true B_i = 0), where a small
  ## denominator inflates CV_i without implying an unstable estimate.
  CV <- sdv / abs(mu)
  CV[!is.finite(CV)] <- Inf
  ## Scale-referenced alternative. Dividing each patient's cross-fold
  ## SD by the TRIAL-level leverage scale rather than by that
  ## patient's own |mean| keeps the statistic dimensionless without
  ## exploding for near-zero-leverage patients. Under arm
  ## cancellation the per-patient CV fires on B_i ~ 0 patients even
  ## though their leverage is estimated just as precisely as everyone
  ## else's, which would disable VOI exactly where it is designed to
  ## help; this variant does not have that failure mode.
  scale_B <- mean(abs(mu), na.rm = TRUE)
  CVs <- if (scale_B > 0) sdv / scale_B else rep(Inf, n)
  list(B = Bmat, CV = CV, CV_scaled = CVs,
       CV90 = as.numeric(quantile(CV, 0.90, na.rm = TRUE, type = 7)),
       CV90_scaled = as.numeric(quantile(CVs, 0.90, na.rm = TRUE, type = 7)),
       frac_near_zero = mean(abs(mu) < 1e-6))
}

## Cross-fold check on the uncertainty score. The supplement's second
## fallback condition compares the input sigma score against a
## cross-fold estimate of the same quantity; here that estimate is the
## per-patient spread over m0 repeat queries in the same-trial
## simulation diagnostic. This is not the frozen primary algorithm.
sigma_score_spearman <- function(m_all, sigma_alloc, m0 = 3, frac_V = 0.30) {
  n <- nrow(m_all)
  m0 <- min(m0, ncol(m_all))
  if (m0 < 2) return(NA_real_)
  nv <- max(round(frac_V * n), 20L)
  V  <- sample.int(n, min(nv, n))
  spread <- apply(m_all[V, seq_len(m0), drop = FALSE], 1, sd)
  if (sd(spread) < .Machine$double.eps || sd(sigma_alloc[V]) < .Machine$double.eps) {
    return(NA_real_)
  }
  suppressWarnings(cor(spread, sigma_alloc[V], method = "spearman"))
}

## Pre-specified fallback rule (Supplement, "Pre-specified threshold"):
##   CV90 <= cv_thresh                      -> proceed with VOI
##   CV90 >  cv_thresh, rho >= rho_thresh   -> fall back to Neyman
##   CV90 >  cv_thresh, rho <  rho_thresh   -> fall back to uniform
voi_diagnostic_decision <- function(CV90, rho_sigma,
                                    cv_thresh = 0.5, rho_thresh = 0.4) {
  if (is.finite(CV90) && CV90 <= cv_thresh) return("VOI")
  if (is.na(rho_sigma) || rho_sigma >= rho_thresh) return("Neyman")
  "Unif"
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
  ## One prespecified partition is shared by the pilot, allocation, and all
  ## final outcome regressions in this replication.
  fold_main <- make_folds(n, L)

  if (identical(sigma_for_alloc, "true")) {
    sigma_alloc <- sigma_LLM_true
  } else {
    sigma_alloc <- sigma_for_alloc
  }

  results <- list()

  ## 1. AIPW
  r <- aipw_xfit(Y, A, X, M = NULL, pi_ = pi_, L = L, include_m = FALSE,
                 fold_id = fold_main)
  results$AIPW <- list(tau = r$tau, se = r$se, K = NA)

  ## 2. CALM (K=1)
  r <- aipw_xfit(Y, A, X, M = m_first, pi_ = pi_, L = L,
                 include_m = TRUE, fold_id = fold_main)
  results$CALM <- list(tau = r$tau, se = r$se, K = rep(1L, n))

  ## 3. Unif (K_i = K)
  K_unif <- rep(as.integer(K_budget), n)
  m_bar  <- avg_predictions(m_all, K_unif)
  r <- aipw_xfit(Y, A, X, M = m_bar, pi_ = pi_, L = L, include_m = TRUE,
                 fold_id = fold_main)
  results$Unif <- list(tau = r$tau, se = r$se, K = K_unif)

  ## 4. Neyman (K_i propto sigma_alloc)
  K_n    <- greedy_allocate_by_fold(sigma_alloc, fold_main, K_budget,
                                    K_min = 1L, K_max = K_max_buffer)
  m_bar  <- avg_predictions(m_all, K_n)
  r <- aipw_xfit(Y, A, X, M = m_bar, pi_ = pi_, L = L, include_m = TRUE,
                 K_alloc = K_n, stratify_by_K = TRUE, fold_id = fold_main)
  results$Neyman <- list(tau = r$tau, se = r$se, K = K_n)

  ## 5. VOI plug-in (K_i propto B_hat * sigma_alloc)
  pilot     <- cross_fit_pilot_B(X, m_first, Y, A, pi_ = pi_, L = L,
                                 fold_id = fold_main, pilot = "separate")
  q_hat_raw <- pilot$B * sigma_alloc

  ## Fit the AIPW estimator for a given allocation score.
  fit_voi <- function(qv) {
    Ka <- greedy_allocate_by_fold(qv, fold_main, K_budget,
                                  K_min = 1L, K_max = K_max_buffer)
    mb <- avg_predictions(m_all, Ka)
    rr <- aipw_xfit(Y, A, X, M = mb, pi_ = pi_, L = L, include_m = TRUE,
                    K_alloc = Ka, stratify_by_K = TRUE,
                    fold_id = fold_main)
    list(tau = rr$tau, se = rr$se, K = Ka)
  }

  ## 5a. VOI as previously reported: winsorized only if winsor_q given
  ##     (default NULL = unstabilized).
  q_hat <- q_hat_raw
  if (!is.null(winsor_q)) {
    q_hat <- winsorize_by_fold(q_hat, fold_main, winsor_q)
  }
  rv <- fit_voi(q_hat)
  results$VOI <- c(rv, list(B_hat = pilot$B))

  ## 5b. VOI oracle (true B_i and true sigma_LLM)
  K_v_or   <- greedy_allocate_by_fold(q_true, fold_main, K_budget,
                                      K_min = 1L, K_max = K_max_buffer)
  m_bar    <- avg_predictions(m_all, K_v_or)
  r <- aipw_xfit(Y, A, X, M = m_bar, pi_ = pi_, L = L, include_m = TRUE,
                 K_alloc = K_v_or, stratify_by_K = TRUE,
                 fold_id = fold_main)
  results$VOI_oracle <- list(tau = r$tau, se = r$se, K = K_v_or,
                             B_true = B_true)

  ## Compute diagnostic and stabilized variants after the six core
  ## estimators so their definitions remain easy to audit.

  ## 5c. VOI stabilized at the prespecified 95th-quantile option.
  q_w  <- winsorize_by_fold(q_hat_raw, fold_main, 0.95)
  r_w  <- fit_voi(q_w)
  results$VOI_w <- r_w

  ## 6. Diagnostic-gated VOI: run the pre-specified pilot-quality
  ##    diagnostic and fall back to Neyman or uniform when it fires.
  ##    The VOI branch uses the stabilized (winsorized) allocation,
  ##    matching the stabilized exploratory analysis.
  disp <- pilot_B_dispersion(
    X, m_first, Y, A, pi_ = pi_, L_diag = 10, pilot = "separate"
  )
  ## The controlled main grid supplies the true synthetic uncertainty
  ## score, so its quality is known by design and no unbudgeted repeat
  ## queries are used for this diagnostic. For externally supplied
  ## scores in a simulation, compare them with the known truth.
  rho_s <- if (identical(sigma_for_alloc, "true")) {
    1
  } else {
    suppressWarnings(cor(sigma_alloc, sigma_LLM_true, method = "spearman"))
  }
  dec   <- voi_diagnostic_decision(disp$CV90, rho_s)
  results$VOI_diag <- switch(dec,
    VOI    = r_w,
    Neyman = results$Neyman,
    Unif   = results$Unif)

  ## Same gate driven by the scale-referenced dispersion instead.
  dec2 <- voi_diagnostic_decision(disp$CV90_scaled, rho_s)
  results$VOI_diag2 <- switch(dec2,
    VOI    = r_w,
    Neyman = results$Neyman,
    Unif   = results$Unif)

  ## 7. VOI-shrink: exploratory stabilization. The shrinkage
  ##    weight comes from the stability diagnostic. Within the
  ##    interpolation lambda = 0 is uniform. The wrapper may instead
  ##    route a zero-weight case to Neyman when the supplied uncertainty
  ##    score is heterogeneous and marked reliable. We record whether
  ##    fallback routed to Neyman or uniform (with constant sigma the
  ##    two coincide anyway).
  lam <- select_lambda(disp$CV90)
  sigma_heterogeneous <- sd(sigma_alloc) > sqrt(.Machine$double.eps)
  fb_route <- if (lam > 0) {
    "none"
  } else if (sigma_heterogeneous && (is.na(rho_s) || rho_s >= 0.4)) {
    "Neyman"
  } else {
    "Unif"
  }
  if (lam == 0 && fb_route == "Neyman") {
    results$VOI_shr <- results$Neyman
  } else {
    K_sh  <- shrinkage_allocate_by_fold(q_hat_raw, fold_main, K_budget, lam,
                                        K_min = 1L, K_max = K_max_buffer)
    m_bar <- avg_predictions(m_all, K_sh)
    r <- aipw_xfit(Y, A, X, M = m_bar, pi_ = pi_, L = L, include_m = TRUE,
                   K_alloc = K_sh, stratify_by_K = TRUE,
                   fold_id = fold_main)
    results$VOI_shr <- list(tau = r$tau, se = r$se, K = K_sh)
  }

  results$.diag <- list(CV90 = disp$CV90,
                        CV90_scaled = disp$CV90_scaled,
                        rho_sigma = rho_s,
                        decision = dec, decision_scaled = dec2,
                        lambda = lam, fb_route = fb_route,
                        frac_near_zero = disp$frac_near_zero)

  attr(results, "X") <- X
  attr(results, "B_true") <- B_true
  attr(results, "sigma_LLM_true") <- sigma_LLM_true
  attr(results, "B_hat") <- pilot$B
  results
}
