## ============================================================
## run_all.R
## Master driver for the CALM-VOI simulations.
##
## Run from the repo root:
##   Rscript scripts/run_all.R          # full run, B = 1000
##   Rscript scripts/run_all.R 100      # quick run, B = 100 (debug)
##   Rscript scripts/run_all.R 500      # B = 500 (still credible MC noise)
##
## Reads no input files. Writes:
##   outputs/results_main.csv          main simulation grid summary
##   outputs/taus_main.rds             raw tau-hat arrays (main grid)
##   outputs/results_sensitivity.csv   sensitivity grid summary
##   outputs/taus_sensitivity.rds      raw tau-hat arrays (sensitivity)
##   outputs/worked_example_mini_mc.csv
##
## Estimated runtime:
##   B = 100:  ~5 min serial
##   B = 500:  ~15 min serial
##   B = 1000: ~30 min serial
## ============================================================

args <- commandArgs(trailingOnly = TRUE)
B <- if (length(args) >= 1) as.integer(args[1]) else 1000L
cat(sprintf("CALM-VOI simulations. B = %d Monte Carlo replications.\n", B))

.find_script <- function(name) {
  for (p in c(paste0("scripts/", name), name)) {
    if (file.exists(p)) return(p)
  }
  stop("Cannot find ", name)
}

cat("\n>>> Step 1/3: main simulation grid (run_main_grid.R)\n\n")
source(.find_script("run_main_grid.R"))
main(B = B)

cat("\n>>> Step 2/3: score-quality sensitivity (run_sensitivity.R)\n\n")
source(.find_script("run_sensitivity.R"))
main_sensitivity(B = B)

cat("\n>>> Step 3/3: worked example (worked_example.R)\n\n")
source(.find_script("worked_example.R"))
main_worked_example()

cat("\nAll done. Output CSVs and RDS files in outputs/\n")
