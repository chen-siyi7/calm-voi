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
## Replications: B = 1000
##
## Outputs:
##   results_main.csv  per-cell summary statistics
##   taus_main.rds     raw tau-hat arrays
##
## Total cells: 6 * 3 * 2 = 36. Estimated runtime: ~15 min serial,
## ~3 min with mclapply on 8 cores (Linux/Mac).
## ============================================================

## Robust path resolution: works from repo root or scripts/
.calm_voi_path <- if (file.exists("R/calm_voi.R")) "R/calm_voi.R" else
                  if (file.exists("../R/calm_voi.R")) "../R/calm_voi.R" else
                  stop("Cannot find R/calm_voi.R. Run from repo root or scripts/.")
source(.calm_voi_path)


if (!dir.exists("outputs")) dir.create("outputs", recursive = TRUE)

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
B_DEFAULT <- 1000
SEED_BASE <- 20260522
ESTIMATORS <- c("AIPW", "CALM", "Unif", "Neyman", "VOI", "VOI_oracle")

## ---- Run one cell ----------------------------------------------
run_cell <- function(dgp_key, dgp_id, sigma_spec, n, K, B, seed_offset) {
  taus <- matrix(NA_real_, nrow = B, ncol = length(ESTIMATORS),
                 dimnames = list(NULL, ESTIMATORS))
  ses  <- matrix(NA_real_, nrow = B, ncol = length(ESTIMATORS),
                 dimnames = list(NULL, ESTIMATORS))
  t0 <- Sys.time()
  for (b in seq_len(B)) {
    set.seed(seed_offset + b)
    r <- one_replication(n = n, K_budget = K, dgp = dgp_id,
                         sigma_eta_spec = sigma_spec)
    for (e in ESTIMATORS) {
      taus[b, e] <- r[[e]]$tau
      ses[b, e]  <- r[[e]]$se
    }
  }
  dt <- as.numeric(Sys.time() - t0, units = "secs")
  row <- data.frame(DGP = dgp_key, n = n, K = K, B = B,
                    runtime_s = round(dt, 1),
                    stringsAsFactors = FALSE)
  for (e in ESTIMATORS) {
    bias <- mean(taus[, e]) - TAU_TRUE
    mc_var <- var(taus[, e])
    mc_se_v <- sqrt(mc_var * (2 / (B - 1)))  ## SE of MC variance (approx)
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
main <- function(B = B_DEFAULT, save_path = "outputs/results_main.csv",
                 taus_path = "outputs/taus_main.rds",
                 verbose = TRUE) {
  cells <- expand.grid(dgp_key = names(DGPS), n = NS, K = KS,
                       stringsAsFactors = FALSE)
  rows <- list()
  all_taus <- list()
  all_ses  <- list()
  for (i in seq_len(nrow(cells))) {
    dgp_key <- cells$dgp_key[i]
    n <- cells$n[i]
    K <- cells$K[i]
    spec <- DGPS[[dgp_key]]
    seed_offset <- SEED_BASE + 1e6 * i
    if (verbose) cat(sprintf("[%d/%d] DGP=%s n=%d K=%d ... ",
                             i, nrow(cells), dgp_key, n, K))
    out <- run_cell(dgp_key, spec$dgp, spec$sigma_spec, n, K, B,
                    seed_offset)
    rows[[i]] <- out$row
    cell_id <- sprintf("%s_n%d_K%d", dgp_key, n, K)
    all_taus[[cell_id]] <- out$taus
    all_ses[[cell_id]]  <- out$ses
    if (verbose) cat(sprintf("runtime %.0fs\n", out$row$runtime_s))
    ## Save partial results to be safe
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
  B <- if (length(args) >= 1) as.integer(args[1]) else B_DEFAULT
  main(B = B)
}
