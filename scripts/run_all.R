## ============================================================
## run_all.R
## Master driver for the CALM-VOI replication package.
##
## Usage:
##   Rscript run_all.R       # recommended submission profile
##   Rscript run_all.R 100   # common-B smoke run
##   CALM_VOI_WORKERS=5 Rscript run_all.R  # parallel grid cells
##
## Recommended profile:
##   main grid and score sensitivity: B = 5000
##   CALM-UQ, external policy, query-noise robustness: B = 2000
##   ACTG query-bank replay: B = 1000
##   turnout query-bank replay: B = 1000
##   robocall query-bank replay: B = 1000
##   IST query-bank replay: B = 200
##   synthetic implementation check: B = 200
##
## Reads the openly licensed IST database below data/ist. Writes:
##   results_main.csv         (main simulation grid summary)
##   taus_main.rds            (raw tau-hat arrays)
##   results_sensitivity.csv  (score-quality sensitivity summary)
##   taus_sensitivity.rds     (raw tau-hat arrays)
##   results_calmuq.csv
##   results_external_pilot.csv
##   results_external_pilot_paired.csv
##   taus_external_pilot.rds
##   results_noise_robustness.csv
##   taus_noise_robustness.rds
##   results_real_application.csv
##   results_real_query_mc.csv
##   results_real_effect_significance.csv
##   results_real_query_comparisons.csv
##   results_turnout_binary_screening.csv
##   results_turnout_binary_allocations.csv
##   results_turnout_binary_replays.csv
##   results_turnout_binary_paired.csv
##   results_robocall_screening.csv
##   results_robocall_allocations.csv
##   results_robocall_replays.csv
##   results_robocall_paired.csv
##   results_ist_application.csv
##   results_ist_query_mc.csv
##   results_ist_effect_significance.csv
##   results_ist_query_comparisons.csv
##   worked_example_mini_mc.csv
##   table_*_rows.tex
##   fig_variance_reduction.pdf
##   fig_allocation_dgp3.pdf
##   fig_sampling_dgp3.pdf
##   fig_sensitivity.pdf
##   fig_real_allocation.pdf
##   fig_ist_allocation.pdf
##
## The recommended serial run can take several hours. A common-B
## argument is intended for smoke tests, not submission results.
## ============================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) {
  stop("Supply at most one common-B override")
}
if (!length(args)) {
  B_main <- 5000L
  B_sensitivity <- 5000L
  B_calmuq <- 2000L
  B_external <- 2000L
  B_noise <- 2000L
  B_query <- 1000L
  B_turnout_query <- 1000L
  B_robocall_query <- 1000L
  B_ist_query <- 200L
  B_pipeline <- 200L
} else {
  B_raw <- suppressWarnings(as.numeric(args[1]))
  if (length(B_raw) != 1L || !is.finite(B_raw) ||
      B_raw < 2 || B_raw != round(B_raw)) {
    stop("The common replication count must be an integer of at least 2")
  }
  B_main <- B_sensitivity <- B_calmuq <- B_external <- B_noise <-
    as.integer(B_raw)
  B_query <- max(20L, min(1000L, as.integer(B_raw)))
  B_turnout_query <- max(20L, min(1000L, as.integer(B_raw)))
  B_robocall_query <- max(20L, min(1000L, as.integer(B_raw)))
  B_ist_query <- max(20L, min(200L, as.integer(B_raw)))
  B_pipeline <- max(20L, min(200L, as.integer(B_raw)))
}
workers_raw <- suppressWarnings(as.numeric(
  Sys.getenv("CALM_VOI_WORKERS", unset = "1")
))
if (length(workers_raw) != 1L || !is.finite(workers_raw) ||
    workers_raw < 1 || workers_raw != round(workers_raw)) {
  stop("CALM_VOI_WORKERS must be a positive integer")
}
workers <- as.integer(workers_raw)
cat(sprintf(
  paste0(
    "CALM-VOI replication profile: main=%d, sensitivity=%d, ",
    "CALM-UQ=%d, external=%d, noise=%d, ACTG replay=%d, ",
    "turnout replay=%d, robocall replay=%d, IST replay=%d, implementation=%d, ",
    "grid workers=%d.\n"
  ),
  B_main, B_sensitivity, B_calmuq, B_external, B_noise,
  B_query, B_turnout_query, B_robocall_query, B_ist_query, B_pipeline, workers
))

cat("\n>>> Step 1/13: main simulation grid (run_main_grid.R)\n\n")
source("run_main_grid.R")
main(B = B_main, workers = workers)

cat("\n>>> Step 2/13: score-quality sensitivity (run_sensitivity.R)\n\n")
source("run_sensitivity.R")
main_sensitivity(B = B_sensitivity, workers = workers)

cat("\n>>> Step 3/13: fixed-gamma CALM-UQ experiment (run_calmuq.R)\n\n")
source("run_calmuq.R")
main_calmuq(B = B_calmuq)

cat("\n>>> Step 4/13: externally trained frozen policy (run_external_pilot.R)\n\n")
source("run_external_pilot.R")
main_external_pilot(B = B_external)

cat("\n>>> Step 5/13: query-noise robustness (run_noise_robustness.R)\n\n")
source("run_noise_robustness.R")
main_noise_robustness(B = B_noise)

cat("\n>>> Step 6/13: worked example (worked_example.R)\n\n")
source("worked_example.R")
main_worked_example(B_mini = B_pipeline)

cat("\n>>> Step 7/13: ACTG 175 application (run_real_application.R)\n\n")
source("run_real_application.R")
main_real_application(B_query = B_query)

cat("\n>>> Step 8/13: household turnout application\n\n")
source("run_turnout_screen.R")

cat("\n>>> Step 9/13: robocall placebo experiment\n\n")
source("run_robocall_screen.R")

cat("\n>>> Step 10/13: International Stroke Trial application\n\n")
source("run_ist_application.R")
main_ist_application(B_query = B_ist_query)

cat("\n>>> Step 11/13: generate figures (make_plots.R)\n\n")
source("make_plots.R")
main_plots()

cat("\n>>> Step 12/13: generate LaTeX rows (make_tables.R)\n\n")
source("make_tables.R")
main_tables()

cat("\n>>> Step 13/13: validate generated outputs (validate_outputs.R)\n\n")
source("validate_outputs.R")
validate_outputs()

cat("\nAll done. Output files in working directory.\n")
