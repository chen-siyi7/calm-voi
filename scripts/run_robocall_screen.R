## Outcome-blind screen of the Yale treatment-versus-placebo robocall trial.

suppressPackageStartupMessages({
  library(foreign)
  library(rpart)
})
source("calm_voi.R")
set.seed(20260814)

raw <- read.dta("data_michigan/D034F02.dta", convert.factors = FALSE)
raw <- raw[
  !is.na(raw$p2008_fi) & raw$hhcount == 1 & raw$treatmen %in% c(1, 2),
]
stopifnot(nrow(raw) == 16805L, sum(raw$treatmen == 1) == 8448L)
as_binary <- function(x) as.numeric(as.character(x))
dat <- data.frame(
  A = as.integer(raw$treatmen == 1),
  Y = as.numeric(raw$p2008_fi),
  absentee_history = raw$av,
  democratic_score = raw$dem,
  area_code = as.numeric(raw$ac),
  house_district = as.numeric(raw$hd),
  primary_2002 = as_binary(raw$p2002),
  primary_2004 = as_binary(raw$p2004),
  general_2002 = as_binary(raw$g2002),
  general_2004 = as_binary(raw$g2004),
  presidential_primary_2008 = as_binary(raw$pp2008)
)
stopifnot(!anyNA(dat))

pi_known <- 0.5
design_fraction <- 0.20
tree_count <- 200L
query_reps <- if (exists("B_robocall_query", inherits = TRUE)) {
  as.integer(get("B_robocall_query", inherits = TRUE))
} else {
  1000L
}
if (!is.finite(query_reps) || query_reps < 20L) {
  stop("The robocall query replay requires at least 20 banks")
}
K_bar <- 5L
K_max <- 10L
baseline_names <- setdiff(names(dat), c("A", "Y"))

design_id <- unlist(lapply(split(seq_len(nrow(dat)), dat$A), function(ii) {
  sample(ii, floor(design_fraction * length(ii)))
}), use.names = FALSE)
analysis_id <- setdiff(seq_len(nrow(dat)), design_id)
design <- dat[design_id, , drop = FALSE]
analysis <- dat[analysis_id, , drop = FALSE]

X_design_raw <- as.matrix(design[, baseline_names, drop = FALSE])
X_analysis_raw <- as.matrix(analysis[, baseline_names, drop = FALSE])
x_center <- colMeans(X_design_raw)
x_scale <- apply(X_design_raw, 2, sd)
x_scale[!is.finite(x_scale) | x_scale == 0] <- 1
X_design <- sweep(sweep(X_design_raw, 2, x_center, "-"), 2, x_scale, "/")
X_analysis <- sweep(sweep(X_analysis_raw, 2, x_center, "-"), 2, x_scale, "/")

pred_design_oob <- matrix(NA_real_, nrow(design), tree_count)
pred_analysis <- matrix(NA_real_, nrow(analysis), tree_count)
mtry <- max(3L, floor(sqrt(ncol(X_design))))
for (b in seq_len(tree_count)) {
  inbag <- sample.int(nrow(design), nrow(design), replace = TRUE)
  oob <- setdiff(seq_len(nrow(design)), unique(inbag))
  vars <- sample(colnames(X_design), mtry, replace = FALSE)
  fit <- rpart(
    Y ~ ., data = data.frame(Y = design$Y[inbag],
                             X_design[inbag, vars, drop = FALSE]),
    method = "anova",
    control = rpart.control(
      minsplit = 80L, minbucket = 30L, cp = 0.0005,
      maxdepth = 7L, xval = 0L
    )
  )
  pred_analysis[, b] <- predict(
    fit, newdata = data.frame(X_analysis[, vars, drop = FALSE])
  )
  if (length(oob)) {
    pred_design_oob[oob, b] <- predict(
      fit, newdata = data.frame(X_design[oob, vars, drop = FALSE])
    )
  }
}
if (any(rowSums(is.finite(pred_design_oob)) < 40L)) {
  stop("Too few OOB tree predictions")
}

m_design <- rowMeans(pred_design_oob, na.rm = TRUE)
m_analysis <- rowMeans(pred_analysis)
sigma <- apply(pred_analysis, 1, sd)
m_center <- mean(m_design)
m_scale <- sd(m_design)
modifier_names <- c(
  "absentee_history", "democratic_score", "primary_2004",
  "presidential_primary_2008"
)
modifier_cols <- match(modifier_names, colnames(X_design))

augmented_matrix <- function(X, m = NULL) {
  if (is.null(m)) return(cbind(`(Intercept)` = 1, X))
  mz <- (m - m_center) / m_scale
  cbind(
    `(Intercept)` = 1, X, m = mz, m2 = mz^2,
    sweep(X[, modifier_cols, drop = FALSE], 1, mz, "*")
  )
}
ridge_fit <- function(Z, y, ridge = 1e-4, maxit = 50L) {
  penalty <- diag(ncol(Z))
  penalty[1, 1] <- 0
  scale <- mean(colSums(Z^2)[-1])
  th <- numeric(ncol(Z))
  th[1] <- qlogis((sum(y) + 0.5) / (length(y) + 1))
  for (iter in seq_len(maxit)) {
    eta <- pmax(-20, pmin(20, as.numeric(Z %*% th)))
    mu <- plogis(eta)
    w <- pmax(mu * (1 - mu), 1e-5)
    working <- eta + (y - mu) / w
    P <- crossprod(Z, w * Z)
    th_new <- drop(solve(
      P + ridge * scale * penalty,
      crossprod(Z, w * working)
    ))
    if (max(abs(th_new - th)) < 1e-8) return(th_new)
    th <- th_new
  }
  th
}
derivative <- function(th, X, m) {
  mz <- (m - m_center) / m_scale
  p <- ncol(X)
  th_m <- th[2 + p]
  th_m2 <- th[3 + p]
  th_int <- th[(4 + p):(3 + p + length(modifier_cols))]
  d_eta <- as.numeric(
    th_m + 2 * th_m2 * mz +
      X[, modifier_cols, drop = FALSE] %*% th_int
  ) / m_scale
  prob <- plogis(as.numeric(augmented_matrix(X, m) %*% th))
  prob * (1 - prob) * d_eta
}

Z_design <- augmented_matrix(X_design, m_design)
theta <- lapply(0:1, function(a) {
  ii <- which(design$A == a)
  ridge_fit(Z_design[ii, , drop = FALSE], design$Y[ii])
})
beta0 <- derivative(theta[[1]], X_analysis, m_analysis)
beta1 <- derivative(theta[[2]], X_analysis, m_analysis)
B <- abs(beta1 + beta0)

n_analysis <- nrow(analysis)
total_budget <- n_analysis * K_bar
K_uniform <- rep(K_bar, n_analysis)
K_uncertainty <- greedy_allocate(
  sigma, total_budget, K_min = 1L, K_max = K_max
)
score_voi_50 <- (0.5 * mean(B) + 0.5 * B) * sigma
K_voi_50 <- greedy_allocate(
  score_voi_50, total_budget, K_min = 1L, K_max = K_max
)
score_voi <- B * sigma
K_voi <- greedy_allocate(
  score_voi, total_budget, K_min = 1L, K_max = K_max
)
stopifnot(sum(K_uncertainty) == total_budget,
          sum(K_voi_50) == total_budget, sum(K_voi) == total_budget)

set.seed(20260815)
fold_id <- integer(n_analysis)
for (a in 0:1) {
  ii <- which(analysis$A == a)
  fold_id[ii] <- sample(rep(seq_len(5L), length.out = length(ii)))
}
aipw <- function(M = NULL) {
  Z <- augmented_matrix(X_analysis, M)
  mu0 <- mu1 <- numeric(n_analysis)
  for (fold in seq_len(5L)) {
    te <- which(fold_id == fold)
    tr <- which(fold_id != fold)
    for (a in 0:1) {
      tr_a <- tr[analysis$A[tr] == a]
      th <- ridge_fit(Z[tr_a, , drop = FALSE], analysis$Y[tr_a])
      pr <- plogis(as.numeric(Z[te, , drop = FALSE] %*% th))
      if (a == 0) mu0[te] <- pr else mu1[te] <- pr
    }
  }
  psi <- analysis$A * (analysis$Y - mu1) / pi_known -
    (1 - analysis$A) * (analysis$Y - mu0) / (1 - pi_known) + mu1 - mu0
  list(tau = mean(psi), se = sd(psi) / sqrt(n_analysis))
}

draw_query_average <- function(K, seed) {
  set.seed(seed)
  max_k <- max(K)
  jj <- matrix(sample.int(tree_count, n_analysis * max_k, replace = TRUE),
               nrow = n_analysis, ncol = max_k)
  running <- out <- numeric(n_analysis)
  for (k in seq_len(max_k)) {
    running <- running + pred_analysis[cbind(seq_len(n_analysis), jj[, k])]
    ii <- which(K == k)
    if (length(ii)) out[ii] <- running[ii] / k
  }
  out
}

allocation_list <- list(
  one_query = rep(1L, n_analysis),
  uniform = K_uniform,
  uncertainty_only = K_uncertainty,
  frozen_voi_50 = K_voi_50,
  frozen_voi = K_voi
)
baseline <- aipw()
full_service <- aipw(m_analysis)
replay <- vector("list", length(allocation_list))
names(replay) <- names(allocation_list)
replay_cores <- if (.Platform$OS.type == "windows") 1L else 4L
for (method in names(allocation_list)) {
  K_method <- allocation_list[[method]]
  rows <- parallel::mclapply(
    seq_len(query_reps),
    function(r) {
      M <- draw_query_average(K_method, 910000L + r)
      fit <- aipw(M)
      c(
        estimate = fit$tau,
        se = fit$se,
        measurement_mse = mean((M - m_analysis)^2),
        squared_error = (fit$tau - full_service$tau)^2
      )
    },
    mc.cores = replay_cores, mc.set.seed = FALSE
  )
  replay[[method]] <- do.call(rbind, rows)
}

Z_analysis <- augmented_matrix(X_analysis, m_analysis)
theta_post <- lapply(0:1, function(a) {
  ii <- which(analysis$A == a)
  ridge_fit(Z_analysis[ii, , drop = FALSE], analysis$Y[ii])
})
beta0_post <- derivative(theta_post[[1]], X_analysis, m_analysis)
beta1_post <- derivative(theta_post[[2]], X_analysis, m_analysis)
B_post <- abs(beta1_post + beta0_post)

cv <- function(x) sd(x) / mean(x)
screening <- data.frame(
  n_total = nrow(dat), n_design = nrow(design), n_analysis = nrow(analysis),
  analysis_treated = sum(analysis$A), query_reps = query_reps,
  baseline_tau = baseline$tau, baseline_se = baseline$se,
  full_service_tau = full_service$tau,
  service_r2_design = cor(m_design, design$Y)^2,
  service_r2_analysis = cor(m_analysis, analysis$Y)^2,
  sigma_cv = cv(sigma), leverage_cv = cv(B), score_cv = cv(B * sigma),
  score_rank_correlation = cor(sigma, B * sigma, method = "spearman"),
  leverage_transport_spearman = cor(B, B_post, method = "spearman"),
  allocation_disagreement = mean(K_uncertainty != K_voi)
)
allocations <- do.call(rbind, lapply(names(allocation_list), function(method) {
  K <- allocation_list[[method]]
  data.frame(
    method = method, K_mean = mean(K), K_sd = sd(K),
    K_p10 = unname(quantile(K, 0.10)),
    K_p50 = unname(quantile(K, 0.50)),
    K_p90 = unname(quantile(K, 0.90)),
    leading_criterion = mean((B * sigma)^2 / K)
  )
}))
replays <- do.call(rbind, lapply(names(replay), function(method) {
  x <- replay[[method]]
  data.frame(
    method = method,
    B_query = query_reps,
    full_service_estimate = full_service$tau,
    mean_estimate = mean(x[, "estimate"]),
    query_bias = mean(x[, "estimate"]) - full_service$tau,
    query_sd = sd(x[, "estimate"]),
    query_sd_mcse = sd(x[, "estimate"]) / sqrt(2 * (query_reps - 1)),
    query_rmse_from_full_service = sqrt(mean(x[, "squared_error"])),
    mean_if_se = mean(x[, "se"]),
    mean_measurement_mse = mean(x[, "measurement_mse"]),
    query_variance_share = var(x[, "estimate"]) / baseline$se^2
  )
}))
paired <- do.call(rbind, lapply(
  c("uniform", "frozen_voi_50", "frozen_voi"),
  function(method) {
    d <- replay[[method]][, "squared_error"] -
      replay[["uncertainty_only"]][, "squared_error"]
    data.frame(
      method = method, mean_squared_error_difference = mean(d),
      paired_mcse = sd(d) / sqrt(length(d)),
      lower_95 = mean(d) - qnorm(0.975) * sd(d) / sqrt(length(d)),
      upper_95 = mean(d) + qnorm(0.975) * sd(d) / sqrt(length(d))
    )
  }
))

write.csv(screening, "results_robocall_screening.csv", row.names = FALSE)
write.csv(allocations, "results_robocall_allocations.csv", row.names = FALSE)
write.csv(replays, "results_robocall_replays.csv", row.names = FALSE)
write.csv(paired, "results_robocall_paired.csv", row.names = FALSE)
saveRDS(
  list(
    screening = screening, allocations = allocations,
    replay_summary = replays, paired = paired, replay = replay,
    split = list(design_id = design_id, analysis_id = analysis_id),
    scores = data.frame(
      sigma = sigma, leverage = B, voi = B * sigma,
      K_uncertainty = K_uncertainty, K_voi = K_voi
    )
  ),
  "results_robocall_screen.rds"
)

print(screening)
print(allocations)
print(replays)
print(paired)
