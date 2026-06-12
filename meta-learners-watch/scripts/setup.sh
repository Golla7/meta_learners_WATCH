#!/bin/bash

# scripts/setup.sh — initialize the project and check dependencies.
# RUN FROM THE REPO ROOT:  ./scripts/setup.sh

echo "=========================================="
echo "META-LEARNERS WATCH — PROJECT SETUP"
echo "=========================================="
echo ""

echo "Creating output directories..."
mkdir -p logs temp results
echo "  logs/    worker log files"
echo "  temp/    per-worker intermediate .rds results"
echo "  results/ final summaries, figures and reports"
echo ""

echo "Setting execution permissions..."
chmod +x scripts/setup.sh scripts/run_parallel.sh scripts/run_analysis.sh 2>/dev/null
echo "  done"
echo ""

echo "Checking required R packages (CRAN)..."
Rscript -e "
packages <- c('dplyr', 'tidyr', 'ggplot2', 'gridExtra', 'viridis', 'scales',
              'forcats', 'ggrepel',                          # plotting
              'SuperLearner', 'nnls', 'caret', 'glmnet',     # nuisance models
              'party', 'partykit', 'model4you', 'ranger',    # forests (Dandl, VI, light config)
              'coin', 'permimp', 'sandwich', 'Hmisc')        # tests, VI, Somers' D

missing <- packages[!packages %in% installed.packages()[, 'Package']]

if (length(missing) > 0) {
  cat('Missing CRAN packages:\n')
  cat(paste(' -', missing, collapse = '\n'), '\n\n')
  cat('Install them with:\n')
  cat('  install.packages(c(', paste0('\"', missing, '\"', collapse = ', '), '))\n')
} else {
  cat('  All required CRAN packages are installed\n')
}
cat('\nOptional: ggridges (extra ridge plot in scripts/aggregate.R)\n')
"

echo ""
echo "Checking GitHub packages..."
Rscript -e "
ok <- TRUE
if (!requireNamespace('benchtm', quietly = TRUE)) {
  ok <- FALSE
  cat('  MISSING: benchtm (Sun et al. 2022 benchmark DGP) — REQUIRED\n')
  cat('    devtools::install_github(\"Sophie-Sun/benchtm\")\n')
}
if (!requireNamespace('rlearner', quietly = TRUE)) {
  ok <- FALSE
  cat('  MISSING: rlearner (R-learner LASSO) — needed for the R-learner method\n')
  cat('    devtools::install_github(\"xnie/rlearner\")\n')
}
if (ok) cat('  benchtm and rlearner are installed\n')
"

echo ""
echo "=========================================="
echo "SETUP COMPLETE"
echo "=========================================="
echo ""
echo "Next steps:"
echo "1. (optional) Adjust parameters in scripts/run_parallel.sh:"
echo "   N_CORES, N_SIMS_TOTAL, SCENARIOS, BETA_LEVELS, METHODS_TO_RUN"
echo "2. Run the simulation (from the repo root):"
echo "   caffeinate ./scripts/run_parallel.sh     # macOS (prevents sleep)"
echo "   ./scripts/run_parallel.sh                # Linux"
echo "3. The controller runs scripts/aggregate.R automatically at the end."
echo "   For the optional deep-dive reports:"
echo "   ./scripts/run_analysis.sh"
echo ""
