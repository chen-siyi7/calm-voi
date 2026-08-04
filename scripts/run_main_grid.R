## ============================================================
## run_main_grid.R
## Main simulation grid for the CALM-VOI manuscript.
##
## DGPs:  1a (sigma=0.5), 1b (sigma=1.0), 1c (sigma=2.0),
##        1d (bimodal: half sigma=0.3, half sigma=2.0),
##        2  (interaction leverage, sigma=1.0),
##        3  (arm cancellation, sigma=1.0)
## Sample sizes: n in {100, 200, 400}
## Budgets:      K in {3, 5}
## Replications: B = 5000
##
## Outputs:
##   results_main.csv  per-cell summary statistics
##   taus_main.rds     raw tau-hat arrays
##
## Total cells: 6 * 3 * 2 = 36. Typical runtime is roughly
## 30--35 minutes serial on a laptop.
## ============================================================

source("calm_voi.R")

DGPS <- list(
  "1a" = list(dgp = "1", sigma_spec = 0.5),
  "1b" = list(dgp = "1", sigma_spec = 1.0),
  "1c" = list(dgp = "1", sigma_spec = 2.0),
  "1d" = list(dgp = "1", sigma_spec = "bimodal"),
  "2"  = list(dgp = "2", sigma_spec = 1.0),
  "3"  = list(dgp = "3", sigma_spec = 1.0)
)
NS <- c(100, 200, 400)
KS <- c(3, 5)
B_DEFAULT <- 5000L
SEED_BASE <- 20260522
## The controlled simulation supplies the true synthetic uncertainty
## scale to isolate allocation-score behavior; it is not an end-to-end
## uncertainty-calibration experiment.
## VOI      = plug-in outcome-leverage score, unstabilized
## VOI_w    = plug-in outcome-leverage score, 95th-quantile winsorized
## VOI_diag = diagnostic-gated, per-patient CV90 rule
## VOI_diag2= diagnostic-gated, scale-referenced CV90 rule
## VOI_shr  = exploratory shrinkage allocation
ESTIMATORS <- c("AIPW", "CALM", "Unif", "Neyman", "VOI", "VOI_w",
                "VOI_oracle", "VOI_diag", "VOI_diag2", "VOI_shr")

## ---- Run one cell ----------------------------------------------
run_cell <- function(dgp_key, dgp_id, sigma_spec, n, K, B, seed_offset) {
  taus <- matrix(NA_real_, nrow = B, ncol = length(ESTIMATORS),
                 dimnames = list(NULL, ESTIMATORS))
  ses  <- matrix(NA_real_, nrow = B, ncol = length(ESTIMATORS),
                 dimnames = list(NULL, ESTIMATORS))
  diag_cv90 <- numeric(B); diag_cv90s <- numeric(B)
  diag_dec  <- character(B); diag_dec2 <- character(B)
  lam_vec   <- numeric(B);   fb_vec    <- character(B)
  t0 <- Sys.time()
  for (b in seq_len(B)) {
    set.seed(seed_offset + b)
    r <- one_replication(n = n, K_budget = K, dgp = dgp_id,
                         sigma_eta_spec = sigma_spec)
    for (e in ESTIMATORS) {
      taus[b, e] <- r[[e]]$tau
      ses[b, e]  <- r[[e]]$se
    }
    diag_cv90[b]  <- r$.diag$CV90
    diag_cv90s[b] <- r$.diag$CV90_scaled
    diag_dec[b]   <- r$.diag$decision
    diag_dec2[b]  <- r$.diag$decision_scaled
    lam_vec[b]    <- r$.diag$lambda
    fb_vec[b]     <- r$.diag$fb_route
  }
  dt <- as.numeric(Sys.time() - t0, units = "secs")
  row <- data.frame(DGP = dgp_key, n = n, K = K, B = B,
                    runtime_s = round(dt, 1),
                    ## diagnostic behavior: median statistic and the
                    ## share of runs each fallback branch was taken
                    diag_CV90_med  = round(median(diag_cv90), 3),
                    diag_CV90s_med = round(median(diag_cv90s), 3),
                    diag_pct_VOI    = round(100 * mean(diag_dec == "VOI"), 1),
                    diag_pct_Neyman = round(100 * mean(diag_dec == "Neyman"), 1),
                    diag_pct_Unif   = round(100 * mean(diag_dec == "Unif"), 1),
                    diag2_pct_VOI   = round(100 * mean(diag_dec2 == "VOI"), 1),
                    ## shrinkage / fallback summary (manuscript table):
                    ## mean shrinkage weight, fallback rate, and split
                    lambda_bar  = round(mean(lam_vec), 3),
                    pr_fb       = round(mean(fb_vec != "none"), 3),
                    pr_fb_unif  = round(mean(fb_vec == "Unif"), 3),
                    pr_fb_ney   = round(mean(fb_vec == "Neyman"), 3),
                    stringsAsFactors = FALSE)
  for (e in ESTIMATORS) {
    bias <- mean(taus[, e]) - TAU_TRUE
    mc_var <- var(taus[, e])
    ## Fourth-moment plug-in MCSE; the normal-theory shortcut can be
    ## seriously anti-conservative in heavy-tailed cells.
    mc_se_v <- mcse_sample_variance(taus[, e])
    lo <- taus[, e] - 1.96 * ses[, e]
    hi <- taus[, e] + 1.96 * ses[, e]
    cov <- mean(lo <= TAU_TRUE & TAU_TRUE <= hi)
    mean_se <- mean(ses[, e])
    row[[paste0(e, "_bias")]]    <- bias
    row[[paste0(e, "_var")]]     <- mc_var
    row[[paste0(e, "_var_mcse")]] <- mc_se_v
    row[[paste0(e, "_cov")]]     <- cov
    row[[paste0(e, "_mean_se")]] <- mean_se
  }
  list(row = row, taus = taus, ses = ses)
}

## ---- Main loop -------------------------------------------------
main <- function(B = B_DEFAULT, save_path = "results_main.csv",
                 taus_path = "taus_main.rds",
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
  cells <- expand.grid(dgp_key = names(DGPS), n = NS, K = KS,
                       stringsAsFactors = FALSE)
  run_one <- function(i) {
    dgp_key <- cells$dgp_key[i]
    n <- cells$n[i]
    K <- cells$K[i]
    spec <- DGPS[[dgp_key]]
    seed_offset <- SEED_BASE + 1e6 * i
    if (verbose) cat(sprintf("[%d/%d] DGP=%s n=%d K=%d ... ",
                             i, nrow(cells), dgp_key, n, K))
    out <- run_cell(dgp_key, spec$dgp, spec$sigma_spec, n, K, B,
                    seed_offset)
    cell_id <- sprintf("%s_n%d_K%d", dgp_key, n, K)
    if (verbose) cat(sprintf("runtime %.0fs\n", out$row$runtime_s))
    list(row = out$row, taus = out$taus, ses = out$ses,
         cell_id = cell_id)
  }
  indices <- seq_len(nrow(cells))
  if (workers > 1L) {
    completed <- parallel::mclapply(
      indices, run_one, mc.cores = workers, mc.preschedule = FALSE
    )
  } else {
    completed <- lapply(indices, run_one)
  }
  rows <- list()
  all_taus <- list()
  all_ses <- list()
  for (i in indices) {
    out <- completed[[i]]
    if (inherits(out, "try-error")) {
      stop(sprintf("simulation cell %d failed: %s", i, out))
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
  B <- if (length(args) >= 1) as.integer(args[1]) else B_DEFAULT
  workers <- if (length(args) >= 2) as.integer(args[2]) else 1L
  main(B = B, workers = workers)
}
