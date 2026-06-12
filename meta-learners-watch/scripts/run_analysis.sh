#!/bin/bash

# scripts/run_analysis.sh — full analysis chain on an existing simulation run.
# RUN FROM THE REPO ROOT, after scripts/run_parallel.sh has finished:
#   ./scripts/run_analysis.sh
#
# 1. aggregate.R             WATCH Objective 1/2/3 summaries + per-scenario PDFs
# 2. analyze_focus_methods.R focused cross-scenario report (results/focus_analysis/)
# 3. plot_paper_figures.R    Sechidis et al. Figure 2/3 reproductions
#                            (results/paper_figures/)

set -e

if [ ! -d temp ] || [ -z "$(ls temp/*.rds 2>/dev/null)" ]; then
    echo "ERROR: no simulation output in temp/. Run ./scripts/run_parallel.sh first."
    exit 1
fi

echo "==> 1/3 aggregate.R (WATCH framework analysis)"
Rscript scripts/aggregate.R

echo ""
echo "==> 2/3 analyze_focus_methods.R (focused cross-scenario report)"
Rscript scripts/analyze_focus_methods.R

echo ""
echo "==> 3/3 plot_paper_figures.R (paper Figure 2/3 reproductions)"
Rscript scripts/plot_paper_figures.R

echo ""
echo "All analysis outputs are in results/."
