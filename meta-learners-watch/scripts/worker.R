#!/usr/bin/env Rscript

# scripts/worker.R
# Parallel simulation worker for the same-footing harness. Run from the REPO
# ROOT (the controller scripts/run_parallel.sh does this automatically).
# Every learner shares one nuisance estimator (5-fold stratified-on-Y
# cross-fit, SL.glmnet+SL.cforest outcomes, propensity estimated even in RCT).
#
# CLI args:
#   worker_id scenario beta start_sim end_sim
#   treatment_scenario methods_to_run_str beta_star_file
#
# Output: temp/results_s{scen}_b{beta}_w{id}_trt-{trt}.rds
# Log:    logs/worker_{id}_s{scen}_b{beta}_trt-{trt}.log

# =============================================================================
# 0. ARGUMENT PARSING
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8) {
  stop(paste0(
    "Usage: scripts/worker.R worker_id scenario beta start_sim end_sim ",
    "treatment_scenario methods_to_run beta_star_file\n",
    "Got ", length(args), " args: ", paste(args, collapse = ", ")
  ))
}

worker_id          <- as.integer(args[1])
scenario           <- as.integer(args[2])
beta_level         <- as.numeric(args[3])   # as.numeric handles decimals
start_sim          <- as.integer(args[4])
end_sim            <- as.integer(args[5])
treatment_scenario <- args[6]
methods_to_run_str <- args[7]
beta_star_file     <- args[8]

methods_to_run <- trimws(unlist(strsplit(methods_to_run_str, ",")))

# =============================================================================
# 1. LOAD HARNESS + DATA GENERATOR (run from the repo root)
# =============================================================================

if (!file.exists("R/meta_learners.R"))
  stop("R/meta_learners.R not found — run scripts/worker.R from the repo root")
COMP_CLEAN_NO_DEMO <- TRUE           # suppress §5 smoke-test block
source("R/meta_learners.R")          # loads all learners + compare_meta_learners

if (!file.exists("R/data_generation.R"))
  stop("R/data_generation.R not found — run scripts/worker.R from the repo root")
source("R/data_generation.R")        # generate_simulation_data() and helpers

# =============================================================================
# 2. LOAD BETA* CALIBRATION
# =============================================================================

beta_star_values <- NULL
if (beta_star_file != "NONE" && file.exists(beta_star_file)) {
  tryCatch({
    beta_star_values <- readRDS(beta_star_file)
    cat(sprintf("Worker %d: loaded beta* from %s\n", worker_id, beta_star_file))
  }, error = function(e) {
    cat(sprintf("WARNING (Worker %d): could not load %s — using defaults\n",
                worker_id, beta_star_file))
  })
} else {
  cat(sprintf("Worker %d: no beta* file — generate_simulation_data defaults apply\n",
              worker_id))
}

# =============================================================================
# 3. LOGGING
# =============================================================================

dir.create("logs", showWarnings = FALSE)
dir.create("temp",  showWarnings = FALSE)

log_file <- sprintf("logs/worker_%d_s%d_b%s_trt-%s.log",
                    worker_id, scenario, as.character(beta_level), treatment_scenario)
sink(log_file, append = FALSE)

cat("========================================\n")
cat(sprintf("CLEAN WORKER %d STARTING\n", worker_id))
cat("========================================\n")
cat(sprintf("Scenario:           %d\n", scenario))
cat(sprintf("Beta Level:         %s\n", as.character(beta_level)))
cat(sprintf("Simulation Range:   %d-%d\n", start_sim, end_sim))
cat(sprintf("Treatment Scenario: %s\n", treatment_scenario))
cat(sprintf("Methods:            %s\n", paste(methods_to_run, collapse = ", ")))
cat(sprintf("Start Time:         %s\n", Sys.time()))
cat("========================================\n\n")

# =============================================================================
# 4. BATCH LOOP
# =============================================================================

tryCatch({

  all_metrics  <- list()
  all_pairwise <- list()
  all_vi       <- list()
  all_tevims   <- list()

  n_batch <- end_sim - start_sim + 1
  cat(sprintf("Worker %d: running %d simulations\n\n", worker_id, n_batch))

  for (s in start_sim:end_sim) {

    cat(sprintf("=== Simulation %d ===\n", s))
    seed_tr <- as.integer(s)
    seed_ev <- seed_tr + 999999L

    # ---- training data ----------------------------------------------------
    dtrain <- tryCatch(
      generate_simulation_data(
        n = 500L, p = 30L,
        scenario          = scenario,
        beta_level        = beta_level,
        treatment_scenario = treatment_scenario,
        seed              = seed_tr,
        beta_star_values  = beta_star_values
      ),
      error = function(e) {
        cat("  ERROR generating train data:", e$message, "\n")
        NULL
      }
    )
    if (is.null(dtrain)) { cat("  Skipping sim", s, "(train data failed)\n\n"); next }

    X        <- dtrain %>% dplyr::select(dplyr::starts_with("X"))
    Y        <- dtrain$Y
    trt      <- dtrain$trt
    tau_true <- dtrain$trt_effect   # true individual CATE

    # ---- independent eval data (for Gao EIF) ------------------------------
    deval <- tryCatch(
      generate_simulation_data(
        n = 500L, p = 30L,
        scenario          = scenario,
        beta_level        = beta_level,
        treatment_scenario = treatment_scenario,
        seed              = seed_ev,
        beta_star_values  = beta_star_values
      ),
      error = function(e) {
        cat("  ERROR generating eval data:", e$message, "\n")
        NULL
      }
    )
    if (is.null(deval)) { cat("  Skipping sim", s, "(eval data failed)\n\n"); next }

    X_eval   <- deval %>% dplyr::select(dplyr::starts_with("X"))
    Y_eval   <- deval$Y
    trt_eval <- deval$trt

    # ---- clean driver: one shared nuisance, all methods -------------------
    t_start <- proc.time()[["elapsed"]]
    res <- tryCatch(
      compare_meta_learners(
        X                  = X,
        Y                  = Y,
        trt                = trt,
        methods            = methods_to_run,
        tau_true           = tau_true,
        X_eval             = X_eval,
        Y_eval             = Y_eval,
        trt_eval           = trt_eval,
        seed               = seed_tr,
        treatment_scenario = treatment_scenario,
        light_config       = FALSE,
        verbose            = TRUE,
        compute_vim        = TRUE
      ),
      error = function(e) {
        cat("  ERROR in compare_meta_learners:", e$message, "\n")
        NULL
      }
    )
    t_elapsed <- proc.time()[["elapsed"]] - t_start

    if (is.null(res)) { cat("  Skipping sim", s, "(driver failed)\n\n"); next }
    cat(sprintf("  Sim %d done: %.1f s total\n", s, t_elapsed))

    # ---- tag rows with simulation metadata and accumulate -----------------
    if (!is.null(res$metrics) && nrow(res$metrics) > 0) {
      m <- res$metrics
      m$Simulation <- s
      m$Scenario   <- scenario
      m$Beta       <- beta_level
      all_metrics[[length(all_metrics) + 1]] <- m
    }

    if (!is.null(res$pairwise_relative_errors) &&
        nrow(res$pairwise_relative_errors) > 0) {
      pw <- res$pairwise_relative_errors
      pw$Simulation         <- s
      pw$Scenario           <- scenario
      pw$Beta               <- beta_level
      pw$Treatment_Scenario <- treatment_scenario
      all_pairwise[[length(all_pairwise) + 1]] <- pw
    }

    if (!is.null(res$results)) {
      vi_entry <- lapply(res$results, function(r) r$vi)
      vi_entry <- vi_entry[!sapply(vi_entry, is.null)]
      if (length(vi_entry) > 0)
        all_vi[[paste0("s", s)]] <- vi_entry
    }

    # TE-VIM TMLE scores (Li 2026 Psi_2/Psi_3) — tag and accumulate.
    if (!is.null(res$tevims_scores) && nrow(res$tevims_scores) > 0) {
      tv <- res$tevims_scores
      tv$Simulation         <- s
      tv$Scenario           <- scenario
      tv$Beta               <- beta_level
      tv$Treatment_Scenario <- treatment_scenario
      all_tevims[[length(all_tevims) + 1]] <- tv
    }

    cat("\n")
  }   # end sim loop

  # ==========================================================================
  # 5. COMBINE AND SAVE
  # ==========================================================================

  combined_metrics  <- if (length(all_metrics)  > 0)
    dplyr::bind_rows(all_metrics) else data.frame()
  combined_pairwise <- if (length(all_pairwise) > 0)
    dplyr::bind_rows(all_pairwise) else NULL
  combined_tevims   <- if (length(all_tevims)   > 0)
    dplyr::bind_rows(all_tevims) else NULL

  output <- list(
    metrics                  = combined_metrics,
    vi_scores                = all_vi,
    pairwise_relative_errors = combined_pairwise,
    tevims_scores            = combined_tevims,  # Li 2026 TE-VIM TMLE (Psi_2/Psi_3)
    ate_estimates            = NULL    # not used by the clean harness
  )

  output_file <- sprintf("temp/results_s%d_b%s_w%d_trt-%s.rds",
                         scenario, as.character(beta_level), worker_id, treatment_scenario)
  saveRDS(output, output_file)

  metrics_rows  <- nrow(combined_metrics)
  pairwise_rows <- if (is.null(combined_pairwise)) 0L else nrow(combined_pairwise)
  vi_sets       <- length(all_vi)
  tevims_rows   <- if (is.null(combined_tevims)) 0L else nrow(combined_tevims)

  cat("========================================\n")
  cat("WORKER COMPLETED\n")
  cat("========================================\n")
  cat(sprintf("Output file:     %s\n", output_file))
  cat(sprintf("Metrics rows:    %d\n", metrics_rows))
  cat(sprintf("Pairwise rows:   %d\n", pairwise_rows))
  cat(sprintf("VI sets:         %d\n", vi_sets))
  cat(sprintf("TE-VIM rows:     %d\n", tevims_rows))
  cat(sprintf("End Time:        %s\n", Sys.time()))
  cat("========================================\n")

}, error = function(e) {
  cat("========================================\n")
  cat(sprintf("FATAL ERROR IN WORKER %d\n", worker_id))
  cat("========================================\n")
  cat(sprintf("Error:  %s\n", e$message))
  cat(sprintf("Time:   %s\n", Sys.time()))
  cat("========================================\n")
  err_file <- sprintf("logs/error_worker_%d_s%d_b%s_trt-%s.rds",
                      worker_id, scenario, as.character(beta_level), treatment_scenario)
  saveRDS(list(error = e$message, traceback = traceback()), err_file)
  quit(status = 1)
})

sink()
