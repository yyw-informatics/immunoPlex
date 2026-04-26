#!/usr/bin/env Rscript
# _run_all_benchmarks.R — Master orchestration for immunoPlex Monte Carlo benchmarks
#
# Staged execution protocol:
#   Stage 1 (dev):  200 reps per cell — for code debugging and smoke testing
#   Stage 2 (pub):  Top up to MCSE target — for publication figures
#
# Usage:
#   Rscript _run_all_benchmarks.R [OPTIONS]
#
# Options:
#   --stage       "dev" (200 reps) or "pub" (MCSE-driven top-up). Default: "dev"
#   --benchmarks  Comma-separated list: "0,1,2,3,4" or "all". Default: "all"
#   --n_cores     Number of parallel cores. Default: 1
#   --cache_dir   Cache directory. Default: inst/simulations/cache
#   --skip_figures Skip figure generation if TRUE. Default: FALSE

# ---- Setup ------------------------------------------------------------------

cat("============================================================\n")
cat("  immunoPlex Benchmark Suite\n")
cat("============================================================\n\n")

# Parse command-line arguments
cli_args <- commandArgs(trailingOnly = TRUE)

parse_cli_arg <- function(flag, default) {
  idx <- which(cli_args == flag)
  if (length(idx) > 0 && idx < length(cli_args)) {
    return(cli_args[idx + 1])
  }
  default
}

STAGE       <- parse_cli_arg("--stage",       "dev")
BENCHMARKS  <- parse_cli_arg("--benchmarks",  "all")
N_CORES     <- as.integer(parse_cli_arg("--n_cores", "1"))
EXTRA_N     <- parse_cli_arg("--extra_n",     NULL)
SKIP_FIGS   <- tolower(parse_cli_arg("--skip_figures", "false")) == "true"

# Locate script directory
script_dir <- {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    dirname(normalizePath(sub("--file=", "", file_arg)))
  } else {
    "inst/simulations"
  }
}

CACHE_DIR <- parse_cli_arg("--cache_dir", file.path(script_dir, "cache"))
if (!dir.exists(CACHE_DIR)) dir.create(CACHE_DIR, recursive = TRUE)

# Stage-specific replication counts
# Dev: 200 reps (fast, for debugging)
# Pub: driven by MCSE targets from the plan
#   - Type I error / coverage: MCSE <= 0.005 -> ~1900 reps
#   - Power: MCSE <= 0.02 -> ~625 reps
#   - Bias ratio: MCSE <= 0.05 -> ~400 reps
#   We use the max (1900) for pub stage to satisfy all targets
STAGE_REPS <- list(
  dev = 200L,
  pub = 2000L
)
n_reps <- STAGE_REPS[[STAGE]]
if (is.null(n_reps)) {
  stop("Unknown stage '", STAGE, "'. Use 'dev' or 'pub'.")
}

# Parse benchmark selection
if (BENCHMARKS == "all") {
  run_benchmarks <- c(0L, 1L, 2L, 3L, 4L, 5L, 6L)
} else {
  run_benchmarks <- as.integer(strsplit(BENCHMARKS, ",")[[1]])
}

cat("Stage:      ", STAGE, " (", n_reps, " reps per cell)\n")
cat("Benchmarks: ", paste(run_benchmarks, collapse = ", "), "\n")
cat("Cores:      ", N_CORES, "\n")
cat("Cache:      ", CACHE_DIR, "\n")
cat("Extra-n:    ", if (!is.null(EXTRA_N)) EXTRA_N else "none (full grid)", "\n")
cat("Figures:    ", if (SKIP_FIGS) "SKIPPED" else "enabled", "\n\n")

# ---- Benchmark Dispatch -----------------------------------------------------

# Track timing and results
run_log <- list()

run_benchmark_script <- function(benchmark_id, script_name) {
  script_path <- file.path(script_dir, script_name)

  if (!file.exists(script_path)) {
    cat("[Benchmark ", benchmark_id, "] Script not found: ", script_path, "\n")
    cat("[Benchmark ", benchmark_id, "] SKIPPED (not yet implemented)\n\n")
    return(list(
      benchmark = benchmark_id,
      status    = "skipped",
      reason    = "script not found",
      elapsed   = 0
    ))
  }

  cat("------------------------------------------------------------\n")
  cat("[Benchmark ", benchmark_id, "] Starting: ", script_name, "\n")
  cat("------------------------------------------------------------\n")

  t0 <- proc.time()

  # Build argument vector for the child script
  child_args <- c(
    "--n_reps",    as.character(n_reps),
    "--n_cores",   as.character(N_CORES),
    "--cache_dir", CACHE_DIR
  )

  status <- tryCatch({
    # Source the script in a child environment so it has access to script_dir
    env <- new.env(parent = globalenv())
    env$script_dir <- script_dir
    # Override commandArgs for the child
    child_cli <- c(
      paste0("--n_reps=", n_reps),   # Parsed differently — we pass as pairs
      child_args
    )
    # Simpler: pass via environment variable
    Sys.setenv(BENCHMARK_N_REPS    = n_reps)
    Sys.setenv(BENCHMARK_N_CORES   = N_CORES)
    Sys.setenv(BENCHMARK_CACHE_DIR = CACHE_DIR)
    if (!is.null(EXTRA_N)) Sys.setenv(BENCHMARK_EXTRA_N = EXTRA_N)

    source(script_path, local = env)
    "success"
  }, error = function(e) {
    cat("[Benchmark ", benchmark_id, "] ERROR: ", conditionMessage(e), "\n")
    "error"
  })

  elapsed <- (proc.time() - t0)["elapsed"]
  cat("[Benchmark ", benchmark_id, "] ", toupper(status),
      " in ", round(elapsed / 60, 1), " min\n\n")

  list(
    benchmark = benchmark_id,
    status    = status,
    elapsed   = elapsed
  )
}


# ---- Execute Benchmarks Sequentially ----------------------------------------

benchmark_scripts <- c(
  "0" = "benchmark_preprocessing.R",
  "1" = "benchmark_censoring.R",
  "2" = "benchmark_mcnemar.R",
  "3" = "benchmark_ancova.R",
  "4" = "benchmark_plsda.R",
  "5" = "benchmark_substitution.R",
  "6" = "benchmark_replicates.R"
)

for (bm in run_benchmarks) {
  bm_key <- as.character(bm)
  if (bm_key %in% names(benchmark_scripts)) {
    result <- run_benchmark_script(bm, benchmark_scripts[[bm_key]])
    run_log[[length(run_log) + 1]] <- result
  } else {
    cat("[Benchmark ", bm, "] Unknown benchmark ID. Skipping.\n\n")
  }
}


# ---- Generate Figures -------------------------------------------------------

if (!SKIP_FIGS) {
  plot_scripts <- c(
    "0" = "plot_B0.R",
    "1" = "plot_B1.R",
    "2" = "plot_B2.R",
    "3" = "plot_B3.R",
    "4" = "plot_B4.R",
    "5" = "plot_B5.R",
    "6" = "plot_B6.R"
  )

  cat("\n------------------------------------------------------------\n")
  cat("  Generating Figures\n")
  cat("------------------------------------------------------------\n\n")

  for (bm in run_benchmarks) {
    bm_key <- as.character(bm)
    if (bm_key %in% names(plot_scripts)) {
      plot_path <- file.path(script_dir, plot_scripts[[bm_key]])
      if (!file.exists(plot_path)) {
        cat("[Figures ", bm, "] Script not found: ", plot_path, "\n")
        next
      }
      cat("[Figures ", bm, "] Running ", plot_scripts[[bm_key]], " ...\n")
      tryCatch({
        env <- new.env(parent = globalenv())
        Sys.setenv(BENCHMARK_CACHE_DIR = CACHE_DIR)
        source(plot_path, local = env)
        cat("[Figures ", bm, "] OK\n")
      }, error = function(e) {
        cat("[Figures ", bm, "] ERROR: ", conditionMessage(e), "\n")
      })
    }
  }
} else {
  cat("\nFigures: SKIPPED (--skip_figures)\n")
}


# ---- Summary Report ---------------------------------------------------------

cat("\n============================================================\n")
cat("  Benchmark Suite Summary\n")
cat("============================================================\n\n")

for (entry in run_log) {
  status_icon <- switch(entry$status,
    success = "OK",
    skipped = "SKIP",
    error   = "FAIL",
    "??"
  )
  cat(sprintf("  Benchmark %d: %-7s  %6.1f min",
              entry$benchmark, status_icon, entry$elapsed / 60))
  if (!is.null(entry$reason)) cat("  (", entry$reason, ")")
  cat("\n")
}

total_time <- sum(sapply(run_log, function(x) x$elapsed))
cat(sprintf("\n  Total time: %.1f min\n", total_time / 60))

# Save run log
log_file <- file.path(CACHE_DIR, paste0("run_log_", STAGE, "_",
                                         format(Sys.time(), "%Y%m%d_%H%M%S"),
                                         ".rds"))
saveRDS(run_log, log_file)
cat("  Run log -> ", log_file, "\n")
cat("\nDone.\n")
