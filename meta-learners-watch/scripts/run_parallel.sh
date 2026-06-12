#!/bin/bash

# scripts/run_parallel.sh
# Parallel simulation controller for the same-footing harness.
# RUN FROM THE REPO ROOT:  ./scripts/run_parallel.sh
# (on macOS, wrap with caffeinate to prevent sleep: caffeinate ./scripts/run_parallel.sh)
#
# Launches scripts/worker.R in parallel for every (scenario, beta) setting;
# per-worker results land in temp/ and logs in logs/. When all workers
# finish, scripts/aggregate.R is run automatically.

# ===================================================================
# MAIN CONFIGURATION
# ===================================================================

# --- Core Settings ---
# N_SIMS_TOTAL must be a multiple of N_CORES.
# Per-sim estimate at n=500 with the full 13-method set: ~6 min for the 11
# shared-nuisance methods + ~1-1.5 min EACH for the two self-contained
# DINA Gao (k2) learners (own per-half SuperLearner stage 1) ≈ 8-10 min/sim.
# RAM budget: ~3 GB per worker → on an 18 GB machine a safe cap is 5 workers.
# Rough wall-clock: 12 settings × 10 sims/worker × ~10 min ≈ 20 h.
N_CORES=5           # parallel workers (≈ 3 GB RAM each)
N_SIMS_TOTAL=50     # 10 sims/worker × 5 workers = 50 sims/setting

# --- Simulation Parameters ---
# Sun et al. (2022) Table 1 scenarios. DGP: benchtm (synthpop X, pre-calibrated
# scen_param). Effect modifiers by scenario:
#   1  step in numeric X11  (prognostic 0.5*I(X1='Y')+X11)
#   2  linear in numeric X14 (prognostic X14 - I(X8='N'))
#   3  subgroup I(X14>0.25 & X1='N')  (prognostic I(X1='N')-0.5*X17)
#   4  subgroup I(X14>0.3 | X4='Y')   (prognostic X11-X14)
# NOTE: scenarios 1-4 do NOT use beta_star_RCT.rds (scen_param is pre-calibrated
# by benchtm::get_b()). p=30, n=500 are fixed inside generate_simulation_data.
SCENARIOS=(1 2 3 4)

# Beta level = paper's b1_rel.  0 = null/homogeneous (Type-I error),
# 1 = beta* (~80% power), 2 = 2*beta* (strong heterogeneity).
BETA_LEVELS=(0 1 2)

# --- Treatment Assignment ---
TREATMENT_SCENARIOS=("RCT")

# --- Methods (must match METHOD_REGISTRY names in R/meta_learners.R) ---
# All 13 harness learners. The worker passes this string to
# compare_meta_learners(methods=...) which validates against METHOD_REGISTRY.
METHODS_TO_RUN="DR-learner (2-Model 2SL),DR-learner (Kennedy),DINA Original,DINA Cross-Fit,DINA Gao (GLM k2),DINA Gao (LASSO k2),DINA Centered,DINA Centered Expanded,Dandl Extension,R-learner (LASSO),TMLE (Li 2026 VTE),WATCH (Residual LASSO),WATCH (Residual Risk)"

# ===================================================================

SIMS_PER_WORKER=$((N_SIMS_TOTAL / N_CORES))
if [ $SIMS_PER_WORKER -eq 0 ]; then SIMS_PER_WORKER=1; fi

TOTAL_JOBS=$((${#SCENARIOS[@]} * ${#BETA_LEVELS[@]} * ${#TREATMENT_SCENARIOS[@]}))
CURRENT_JOB=0

echo "============================================"
echo "PARALLEL SIMULATION CONTROLLER"
echo "============================================"
echo "Total sims per setting: $N_SIMS_TOTAL"
echo "Workers:                $N_CORES ($SIMS_PER_WORKER sims/worker)"
echo "Scenarios:              ${SCENARIOS[@]}"
echo "Beta levels:            ${BETA_LEVELS[@]}"
echo "Treatment scenarios:    ${TREATMENT_SCENARIOS[@]}"
echo "Methods:                $METHODS_TO_RUN"
echo "Total jobs:             $TOTAL_JOBS"
echo "Output:                 temp/ and logs/"
echo "============================================"
echo ""

# ===================================================================
# CALIBRATION FILE VALIDATION
# ===================================================================
# data/beta_star_RCT.rds is loaded by scripts/worker.R and forwarded to
# generate_simulation_data() so that legacy scenarios 5-8 work if ever added.
# For the benchtm scenarios 1-4 it is ignored internally (b0/b1 come
# pre-calibrated from benchtm::scen_param), but we validate it exists as a
# safety guard.

echo "============================================"
echo "VERIFYING CALIBRATION FILES"
echo "============================================"

CALIBRATION_MISSING=0
for trt_scenario in "${TREATMENT_SCENARIOS[@]}"; do
    BETA_STAR_FILE="data/beta_star_${trt_scenario}.rds"
    if [ -f "$BETA_STAR_FILE" ]; then
        echo "  Found: $BETA_STAR_FILE"
    else
        echo "  MISSING: $BETA_STAR_FILE"
        CALIBRATION_MISSING=1
    fi
done

if [ $CALIBRATION_MISSING -eq 1 ]; then
    echo ""
    echo "ERROR: Missing calibration file(s)."
    echo "Generate with (in R from repo root):"
    echo "  source('R/data_generation.R')"
    echo "  r <- calculate_beta_star_for_scenarios(scenarios_to_calc=1:4, treatment_scenario='RCT')"
    echo "  saveRDS(r, 'data/beta_star_RCT.rds')"
    exit 1
fi

echo "  All calibration files verified."
echo "============================================"
echo ""

# ===================================================================
# PREPARE OUTPUT DIRS
# ===================================================================
mkdir -p temp logs
# Clear the previous run's outputs
rm -f temp/*.rds
rm -f logs/*.log

# ===================================================================
# MAIN LOOP
# ===================================================================

for trt_scenario in "${TREATMENT_SCENARIOS[@]}"; do
  for scenario in "${SCENARIOS[@]}"; do
    for beta in "${BETA_LEVELS[@]}"; do

      CURRENT_JOB=$((CURRENT_JOB + 1))
      BETA_STAR_FILE="data/beta_star_${trt_scenario}.rds"

      echo "============================================"
      echo "JOB $CURRENT_JOB / $TOTAL_JOBS"
      echo "  Treatment Scenario: $trt_scenario"
      echo "  Data Scenario:      $scenario"
      echo "  Beta Level:         $beta"
      echo "  Start Time:         $(date)"
      echo "============================================"

      # Launch workers in parallel
      for worker_id in $(seq 1 $N_CORES); do
          start_sim=$(((worker_id - 1) * SIMS_PER_WORKER + 1))
          end_sim=$((worker_id * SIMS_PER_WORKER))
          if [ $worker_id -eq $N_CORES ]; then end_sim=$N_SIMS_TOTAL; fi
          if [ $start_sim -gt $end_sim ]; then continue; fi

          echo "  Worker $worker_id: sims $start_sim-$end_sim"

          Rscript scripts/worker.R \
              $worker_id $scenario $beta $start_sim $end_sim \
              "$trt_scenario" "$METHODS_TO_RUN" "$BETA_STAR_FILE" &
      done

      # Wait for all workers for this setting
      wait
      echo "  JOB $CURRENT_JOB DONE at $(date)"
      echo ""

    done
  done
done

echo "============================================"
echo "ALL SIMULATIONS COMPLETE"
echo "============================================"
echo "Running scripts/aggregate.R..."
Rscript scripts/aggregate.R
echo "Done. All outputs are in results/ (summary CSVs + WATCH analysis PDFs)."
