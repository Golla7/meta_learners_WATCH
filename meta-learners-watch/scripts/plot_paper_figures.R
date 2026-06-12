#!/usr/bin/env Rscript
# ============================================================================
# scripts/plot_paper_figures.R  (run from the repo root, after a simulation run)
# Reproduce figures from Sechidis et al. (arXiv:2502.00713v2) that can be built
# from the simulation output already in temp/*.rds — no re-run needed.
#
#   Figure 2  Objective 1(i): ECDF of NULL (beta1=0) p-values. Uniform => the
#             ECDF tracks the diagonal (good type-I behaviour).
#   Figure 3  Objective 1(ii): average p-value vs beta1/beta1*; lower = more
#             powerful test.
#
# Two flavours of each:
#   (A) PAPER-FAITHFUL: one DR-learner, DR-Maximum vs DR-Quadratic statistic
#       (exactly the paper's Fig 2 / Fig 3).
#   (B) ALL-METHODS:    every meta-learner from the run, using the max-type
#       statistic, so you can compare Objective-1 behaviour across methods.
#
# Uses the global-test p-values already stored per simulation (Het_P_Value_Max /
# Het_P_Value_Quad). Fully SOFT-CODED: scenarios, beta levels, labels, excluded
# methods and colours all come from the data + R/scenario_meta.R.
#
# Note: our run uses beta1/beta1* in {0,1,2} (paper sweeps {0,0.5,1,1.5,2}), so
# the Figure-3 lines show 3 points instead of 5.
# ============================================================================
suppressWarnings(suppressMessages({
  library(dplyr); library(tidyr); library(ggplot2)
}))
source("R/scenario_meta.R")   # EXCLUDE_METHODS, METHOD_PAL, scenario_label

# --- CONFIG ----------------------------------------------------------------
# Paper-faithful panels use this single method = the paper's cross-fit
# DR-learner (SuperLearner glmnet+cforest) = "DR-learner (2-Model 2SL)".
PAPER_FIG_METHOD <- "DR-learner (2-Model 2SL)"
OUT <- "results/paper_figures"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
STAT_COLS <- c("DR-Maximum Statistic" = "#D7301F", "DR-Quadratic Statistic" = "#2C7FB8")

# --- pool per-sim metrics ---------------------------------------------------
files <- list.files("temp", pattern = "rds$", full.names = TRUE)
M <- do.call(rbind, lapply(files, function(f) tryCatch(readRDS(f)$metrics, error = function(e) NULL)))
M <- M[!is.na(M$Het_P_Value_Max), , drop = FALSE]

# meta-learners present (those in the shared palette), minus excluded ones.
methods_present <- setdiff(intersect(unique(M$Method), names(METHOD_PAL)), EXCLUDE_METHODS)
M <- M[M$Method %in% methods_present, , drop = FALSE]

scen_levels <- sort(unique(M$Scenario))
scen_lab_lvls <- vapply(scen_levels, scenario_label, character(1))
mk_scenlab <- function(s) factor(vapply(s, scenario_label, character(1)), levels = scen_lab_lvls)
M$ScenLab <- mk_scenlab(M$Scenario)
betas <- sort(unique(M$Beta))

# ============================================================================
# (A) PAPER-FAITHFUL: DR-learner, Max vs Quadratic
# ============================================================================
mm <- M[M$Method == PAPER_FIG_METHOD & !is.na(M$Het_P_Value_Quad), , drop = FALSE]
if (nrow(mm) > 0) {
  n_per <- nrow(mm[mm$Beta == 0, ]) / length(scen_levels)

  # Fig 2 (faithful): null p-value ECDF, Max vs Quad
  f2 <- mm %>% filter(Beta == 0) %>%
    select(ScenLab, `DR-Maximum Statistic` = Het_P_Value_Max,
           `DR-Quadratic Statistic` = Het_P_Value_Quad) %>%
    pivot_longer(-ScenLab, names_to = "Statistic", values_to = "p")
  fig2 <- ggplot(f2, aes(p, color = Statistic)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    stat_ecdf(linewidth = 0.9) + facet_wrap(~ ScenLab, nrow = 1) +
    scale_color_manual(values = STAT_COLS, name = NULL) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(title = "Comparison of the two heterogeneity tests (Sec. 2.2.2) - Objective 1(i)",
         subtitle = sprintf("ECDF of null (beta1=0) p-values; uniform => diagonal.   [%s, %d sims/scenario]",
                            PAPER_FIG_METHOD, round(n_per)),
         x = "p-values", y = "ECDF") +
    theme_bw(base_size = 11) + theme(legend.position = "top")
  ggsave(file.path(OUT, "Figure2_null_pvalue_ECDF.png"), fig2, width = 11, height = 3.3, dpi = 150)

  # Fig 3 (faithful): avg p-value vs beta, Max vs Quad
  f3 <- mm %>% group_by(ScenLab, Beta) %>%
    summarise(`DR-Maximum Statistic` = mean(Het_P_Value_Max, na.rm = TRUE),
              `DR-Quadratic Statistic` = mean(Het_P_Value_Quad, na.rm = TRUE), .groups = "drop") %>%
    pivot_longer(c(`DR-Maximum Statistic`, `DR-Quadratic Statistic`),
                 names_to = "Statistic", values_to = "avg_p")
  fig3 <- ggplot(f3, aes(Beta, avg_p, color = Statistic)) +
    geom_line(linewidth = 0.9) + geom_point(size = 2.4) + facet_wrap(~ ScenLab, nrow = 1) +
    scale_color_manual(values = STAT_COLS, name = NULL) + scale_x_continuous(breaks = betas) +
    labs(title = "Comparison of the two heterogeneity tests (Sec. 2.2.2) - Objective 1(ii)",
         subtitle = sprintf("Average p-value vs heterogeneity strength; lower = more powerful.   [%s]", PAPER_FIG_METHOD),
         x = expression(beta[1] / beta[1]^"*"), y = "Average p-value") +
    theme_bw(base_size = 11) + theme(legend.position = "top")
  ggsave(file.path(OUT, "Figure3_avg_pvalue_vs_beta.png"), fig3, width = 11, height = 3.3, dpi = 150)
}

# ============================================================================
# (B) ALL-METHODS: max-type statistic across every mega-run meta-learner
# ============================================================================
meth_lvls <- intersect(names(METHOD_PAL), methods_present)  # palette order
M$MethodF <- factor(M$Method, levels = meth_lvls)

# Fig 2 (all methods): null p-value ECDF
fig2a <- ggplot(M[M$Beta == 0, ], aes(Het_P_Value_Max, color = MethodF)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  stat_ecdf(linewidth = 0.8) + facet_wrap(~ ScenLab, nrow = 1) +
  scale_color_manual(values = METHOD_PAL, name = "Method", drop = TRUE) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(title = "Objective 1(i): null p-value ECDF — all meta-learners (max-type test)",
       subtitle = "Under beta1=0 a well-calibrated test => ECDF on the diagonal.",
       x = "p-values", y = "ECDF") +
  theme_bw(base_size = 11) + theme(legend.position = "right")
ggsave(file.path(OUT, "Figure2_allmethods_null_ECDF.png"), fig2a, width = 12, height = 3.4, dpi = 150)

# Fig 3 (all methods): avg p-value vs beta
f3a <- M %>% group_by(ScenLab, Beta, MethodF) %>%
  summarise(avg_p = mean(Het_P_Value_Max, na.rm = TRUE), .groups = "drop")
fig3a <- ggplot(f3a, aes(Beta, avg_p, color = MethodF)) +
  geom_line(linewidth = 0.8) + geom_point(size = 2) + facet_wrap(~ ScenLab, nrow = 1) +
  scale_color_manual(values = METHOD_PAL, name = "Method", drop = TRUE) +
  scale_x_continuous(breaks = betas) +
  labs(title = "Objective 1(ii): average p-value vs heterogeneity strength — all meta-learners (max-type)",
       subtitle = "Lower average p-value = more powerful global heterogeneity test.",
       x = expression(beta[1] / beta[1]^"*"), y = "Average p-value") +
  theme_bw(base_size = 11) + theme(legend.position = "right")
ggsave(file.path(OUT, "Figure3_allmethods_avg_pvalue.png"), fig3a, width = 12, height = 3.4, dpi = 150)

# combined PDF (faithful + all-methods)
pdf(file.path(OUT, "Paper_Figures_Objective1.pdf"), width = 12, height = 3.5, onefile = TRUE)
if (exists("fig2")) print(fig2); if (exists("fig3")) print(fig3)
print(fig2a); print(fig3a); invisible(dev.off())

cat(sprintf("Saved Objective-1 figures to %s/\n  methods: %s\n  scenarios: %s | betas: %s\n",
            OUT, paste(methods_present, collapse = ", "),
            paste(scen_levels, collapse = ","), paste(betas, collapse = ",")))
