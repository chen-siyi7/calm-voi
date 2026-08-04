## Reproducible fixed-gamma CALM-UQ experiment for Supplement S9.
## This is deliberately separate from the archived exploratory,
## data-adaptive Python experiment, which uses a different DGP.

source("calm_voi.R")

calmuq_one <- function(gamma, X, A, Y, m_first, u, fold_id, pi_ = 0.5) {
  n <- length(Y)
  mu0 <- mu1 <- numeric(n)
  Z <- cbind(1, X, m_first)
  w <- (1 - u)^gamma
  for (ell in sort(unique(fold_id))) {
    te <- which(fold_id == ell)
    tr <- which(fold_id != ell)
    for (a in 0:1) {
      ii <- tr[A[tr] == a]
      fit <- lm.wfit(Z[ii, , drop = FALSE], Y[ii], w[ii])
      cf <- fit$coefficients
      cf[is.na(cf)] <- 0
      pred <- as.numeric(Z[te, , drop = FALSE] %*% cf)
      if (a == 1) mu1[te] <- pred else mu0[te] <- pred
    }
  }
  score <- A * (Y - mu1) / pi_ -
    (1 - A) * (Y - mu0) / (1 - pi_) + mu1 - mu0
  mean(score)
}

CALMUQ_B_DEFAULT <- 2000L

main_calmuq <- function(B = CALMUQ_B_DEFAULT,
                        save_path = "results_calmuq.csv",
                        verbose = TRUE) {
  if (length(B) != 1L || !is.finite(B) || B < 2 || B != round(B)) {
    stop("B must be a finite integer of at least 2")
  }
  B <- as.integer(B)
  gamma_grid <- c(0, 0.5, 1, 2)
  tau <- matrix(NA_real_, B, length(gamma_grid),
                dimnames = list(NULL, as.character(gamma_grid)))
  for (b in seq_len(B)) {
    set.seed(6200000L + b)
    base <- generate_baseline(400)
    sigma_eta <- make_sigma_eta("bimodal", 400)
    out <- make_outcomes("1", base$X, base$Tstar, base$A)
    m_first <- llm_queries(base$X, base$Tstar, sigma_eta, 1)[, 1]
    fold_id <- make_folds(400, 5)
    ## Map the two oracle uncertainty levels into [0, 0.9], keeping
    ## strictly positive WLS weights at every gamma.
    u <- 0.9 * (sigma_eta - min(sigma_eta)) /
      max(diff(range(sigma_eta)), .Machine$double.eps)
    for (j in seq_along(gamma_grid)) {
      tau[b, j] <- calmuq_one(gamma_grid[j], base$X, base$A, out$Y,
                              m_first, u, fold_id)
    }
    if (verbose && b %% 100L == 0L) cat("CALM-UQ", b, "/", B, "\n")
  }
  v <- apply(tau, 2, var)
  out <- data.frame(
    gamma = gamma_grid,
    B = B,
    bias = colMeans(tau) - TAU_TRUE,
    variance = v,
    variance_mcse = apply(tau, 2, mcse_sample_variance),
    delta_pct = 100 * (v / v[1] - 1)
  )
  write.csv(out, save_path, row.names = FALSE)
  saveRDS(tau, "taus_calmuq.rds")
  invisible(out)
}

if (!interactive() && sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  B <- if (length(args)) as.integer(args[1]) else CALMUQ_B_DEFAULT
  main_calmuq(B)
}
