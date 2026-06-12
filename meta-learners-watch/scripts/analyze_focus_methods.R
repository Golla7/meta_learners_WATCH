#!/usr/bin/env Rscript
# =====================================================================
# scripts/analyze_focus_methods.R
# Focused cross-scenario analysis of the selected HTE methods.
# RUN FROM THE REPO ROOT, AFTER scripts/aggregate.R:
#   Rscript scripts/analyze_focus_methods.R
# Pools per-simulation metrics from temp/*.rds and combines them with the
# results/final_*.csv summaries. Outputs go to results/focus_analysis/.
# =====================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
  library(scales); library(forcats); library(viridis); library(ggrepel)
  library(grid); library(gridExtra)
})

OUT <- "results/focus_analysis"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# Soft-coded scenario metadata + shared analysis config (true modifiers, labels,
# EXCLUDE_METHODS, METHOD_PAL) so the focus analysis never drifts. See R/scenario_meta.R.
source("R/scenario_meta.R")

# ---- focus methods, fixed order & colorblind-friendly palette ----
FOCUS_ALL <- c(
  "DINA Original", "DINA Cross-Fit", "DINA Gao (LASSO k2)",
  "DINA Centered", "DINA Centered Expanded",   # §3z centered-design variants
  "DR-learner (2-Model 2SL)",  "R-learner (LASSO)",
  "Dandl Extension"
)
FOCUS <- setdiff(FOCUS_ALL, EXCLUDE_METHODS)   # EXCLUDE_METHODS from R/scenario_meta.R
PAL   <- METHOD_PAL                            # shared palette (R/scenario_meta.R)
# CATE-error methods that participate in pairwise comparisons (VTE test excluded)
FOCUS_CATE <- setdiff(FOCUS, "TMLE (Li 2026 VTE)")

# Soft scenario labels derived from the DGP (R/scenario_meta.R::scenario_label),
# e.g. "S1 (X11)". Computed once per unique scenario for speed.
lab_scen <- function(s) {
  s  <- as.numeric(s)
  us <- sort(unique(s))
  labs <- setNames(vapply(us, scenario_label, character(1)), as.character(us))
  factor(labs[as.character(s)], levels = labs)
}
mk_method <- function(m) factor(m, levels = FOCUS)

theme_set(theme_bw(base_size = 12))
GUIDE <- guides(color = guide_legend(ncol = 1), fill = guide_legend(ncol = 1))

# ---------------------------------------------------------------------
# 1. Pool per-simulation metrics from temp/*.rds
# ---------------------------------------------------------------------
cat("Pooling per-simulation metrics from temp/ ...\n")
files <- list.files("temp", pattern = "\\.rds$", full.names = TRUE)
raw <- lapply(files, function(f) {
  m <- tryCatch(readRDS(f)$metrics, error = function(e) NULL)
  if (is.null(m)) return(NULL)
  m[m$Method %in% FOCUS, , drop = FALSE]
})
raw <- bind_rows(raw)
cat(sprintf("  pooled %d rows from %d files\n", nrow(raw), length(files)))

# Clean harness stores Het_Detected_Max / Het_Detected_Quad (no plain
# Het_Detected). Alias the max-type detection to the legacy column name so the
# downstream summarise() works unchanged.
if (!"Het_Detected" %in% names(raw) && "Het_Detected_Max" %in% names(raw))
  raw$Het_Detected <- raw$Het_Detected_Max

raw <- raw %>%
  mutate(Method = mk_method(Method), Scen = lab_scen(Scenario))

# de-duplicate sims (worker batches share Simulation index): give each a unique id
raw <- raw %>% group_by(Scenario, Beta, Method) %>%
  mutate(sim_id = row_number()) %>% ungroup()

n_per_cell <- raw %>% count(Scenario, Beta, Method) %>% summarise(rng = paste(range(n), collapse = "-")) %>% pull(rng)
cat("  sims per (scenario,beta,method):", n_per_cell, "\n")

# ---------------------------------------------------------------------
# 2. Per-cell summary with Monte-Carlo SE (mean +/- 1.96*SD/sqrt(n))
# ---------------------------------------------------------------------
se <- function(x) { x <- x[is.finite(x)]; if (length(x) < 2) return(NA_real_); sd(x)/sqrt(length(x)) }
mn <- function(x) mean(x[is.finite(x)], na.rm = TRUE)

cell <- raw %>%
  group_by(Scenario, Scen, Beta, Method) %>%
  summarise(
    n          = sum(is.finite(MSE)),
    # NOTE: compute every *_se BEFORE the matching mean — dplyr evaluates
    # summarise() expressions sequentially, so `MSE = mn(MSE)` would rebind
    # the name `MSE` to a scalar and break a later `se(MSE)`.
    MSE_se = se(MSE), MSE = mn(MSE),
    Bias_se = se(Bias), Bias = mn(Bias),
    Var_se = se(Variance), Variance = mn(Variance),
    Cor_se = se(Correlation), Correlation = mn(Correlation),
    Som_se = se(Somers_D), Somers_D = mn(Somers_D),
    Detection = 100 * mn(Het_Detected),
    Detection_Max = 100 * mn(Het_Detected_Max),
    Detection_Quad = 100 * mn(Het_Detected_Quad),
    Time = mn(Time_Sec),
    .groups = "drop"
  )
write.csv(cell, file.path(OUT, "focus_summary_by_cell.csv"), row.names = FALSE)

# ---------------------------------------------------------------------
# helper: line plot with ribbon CI, faceted by scenario
# ---------------------------------------------------------------------
line_facets <- function(df, y, ylo, yhi, ytitle, subtitle, logy = FALSE, hline = NULL) {
  p <- ggplot(df, aes(Beta, .data[[y]], color = Method, fill = Method, group = Method))
  if (!is.null(hline)) p <- p + geom_hline(yintercept = hline, linetype = "dashed", color = "grey40")
  if (all(c(ylo, yhi) %in% names(df)))
    p <- p + geom_ribbon(aes(ymin = .data[[ylo]], ymax = .data[[yhi]]),
                         alpha = 0.12, color = NA)
  p <- p +
    geom_line(linewidth = 0.9) + geom_point(size = 1.9) +
    facet_wrap(~ Scen, scales = "free_y", ncol = 3) +
    scale_color_manual(values = PAL, drop = FALSE) +
    scale_fill_manual(values = PAL, drop = FALSE) +
    labs(x = expression(beta~"(heterogeneity strength)"), y = ytitle,
         title = ytitle, subtitle = subtitle) +
    GUIDE
  if (logy) p <- p + scale_y_log10(labels = label_number())
  p
}

addCI <- function(df, m, s) df %>% mutate(lo = .data[[m]] - 1.96*.data[[s]],
                                          hi = .data[[m]] + 1.96*.data[[s]])

# ---------------------------------------------------------------------
# 3. Build plots
# ---------------------------------------------------------------------
plots <- list()

# (1) MSE
d <- addCI(cell, "MSE", "MSE_se")
plots$MSE <- line_facets(d, "MSE", "lo", "hi",
  "CATE estimation error (MSE)",
  "Lower is better. Ribbons = 95% Monte-Carlo CI. beta=0 is the null (constant effect).")

# (2) |Bias|
d <- addCI(cell, "Bias", "Bias_se")
plots$Bias <- line_facets(d, "Bias", "lo", "hi",
  "CATE bias (signed)",
  "Closer to 0 is better. Dashed line = 0.", hline = 0)

# (3) Variance
d <- addCI(cell, "Variance", "Var_se")
plots$Variance <- line_facets(d, "Variance", "lo", "hi",
  "CATE estimate variance",
  "Lower = more stable across the sample. Ribbons = 95% MC CI.")

# (4) Correlation with true CATE (beta>0 only; undefined at null)
d <- addCI(filter(cell, Beta > 0), "Correlation", "Cor_se")
plots$Correlation <- line_facets(d, "Correlation", "lo", "hi",
  "Correlation with true CATE",
  "Higher is better (ranking quality). Undefined at beta=0, so shown for beta>0.")

# (5) Somers' D
d <- addCI(filter(cell, Beta > 0), "Somers_D", "Som_se")
plots$Somers <- line_facets(d, "Somers_D", "lo", "hi",
  "Somers' D (concordance with truth)",
  "Higher is better. Shown for beta>0.")

# (6) Detection rate (power; at beta=0 it is the type-I error rate)
plots$Detection <- line_facets(cell, "Detection", NA, NA,
  "Heterogeneity detection rate (%)",
  "beta=0 = type-I error (target ~5%); beta>0 = power. Dashed line = 5%.",
  hline = 5)

# (6b) three detection variants for the focus methods, faceted by scenario at beta=2
det_long <- cell %>%
  select(Scen, Beta, Method, Linear = Detection, Max = Detection_Max, Quad = Detection_Quad) %>%
  pivot_longer(c(Linear, Max, Quad), names_to = "Test", values_to = "Rate")
plots$DetectionVariants <- ggplot(filter(det_long, Beta == 2),
    aes(Method, Rate, fill = Test)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  facet_wrap(~ Scen, ncol = 3) +
  scale_fill_brewer(palette = "Set2", name = "Test statistic") +
  coord_flip() +
  labs(title = "Detection rate by test statistic (beta = 2)",
       subtitle = "Linear vs Max vs Quadratic heterogeneity test, strongest signal.",
       x = NULL, y = "Detection rate (%)")

# (7) MSE distribution at strongest signal (beta=2)
plots$MSE_dist <- ggplot(filter(raw, Beta == 2),
    aes(fct_rev(Method), MSE, fill = Method)) +
  geom_boxplot(alpha = 0.85, outlier.size = 0.7) +
  facet_wrap(~ Scen, scales = "free_x", ncol = 3) +
  scale_fill_manual(values = PAL, drop = FALSE) +
  coord_flip() +
  labs(title = "Per-simulation MSE distribution (beta = 2)",
       subtitle = "Spread across the Monte-Carlo replications.",
       x = NULL, y = "MSE") +
  theme(legend.position = "none")

# (8) Cross-scenario heatmap: fill by within-column rank (1 = best), label = mean MSE.
#     Absolute MSE spans 0.003-0.5, so colouring by rank keeps the ordering
#     legible in every cell while the annotation preserves the magnitude.
heat <- cell %>%
  group_by(Scen, Beta) %>%
  mutate(rank = rank(MSE, ties.method = "min")) %>% ungroup()
plots$Heatmap <- ggplot(heat, aes(Scen, fct_rev(Method), fill = rank)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.3f", MSE)), size = 2.5) +
  facet_wrap(~ paste0("beta = ", Beta), ncol = 4) +
  scale_fill_gradientn(colours = c("#1A9850", "#FFFFBF", "#D73027"),
                       name = "MSE rank\n(1 = best)", breaks = c(1, 5, 10)) +
  labs(title = "Cross-scenario MSE: rank within each scenario x beta cell",
       subtitle = "Green = best (lowest MSE), red = worst. Numbers are the mean MSE.",
       x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

# (9) Overall ranking: average rank across all 20 cells (by MSE) + by Correlation
rank_mse <- cell %>% group_by(Scenario, Beta) %>%
  mutate(r = rank(MSE, ties.method = "average")) %>% ungroup() %>%
  group_by(Method) %>% summarise(AvgRank = mean(r), .groups = "drop") %>%
  mutate(Metric = "MSE (lower better)")
rank_cor <- cell %>% filter(Beta > 0) %>% group_by(Scenario, Beta) %>%
  mutate(r = rank(-Correlation, ties.method = "average")) %>% ungroup() %>%
  group_by(Method) %>% summarise(AvgRank = mean(r), .groups = "drop") %>%
  mutate(Metric = "Correlation (higher better)")
rank_df <- bind_rows(rank_mse, rank_cor)
write.csv(rank_df, file.path(OUT, "focus_overall_ranking.csv"), row.names = FALSE)

# Reframe as a 2D tradeoff: magnitude accuracy (MSE rank) vs ranking quality
# (correlation rank). Axes reversed so "better" is up & to the right.
rank_wide <- rank_df %>%
  mutate(Metric = recode(Metric, "MSE (lower better)" = "MSE_rank",
                                  "Correlation (higher better)" = "Cor_rank")) %>%
  pivot_wider(names_from = Metric, values_from = AvgRank)
plots$Ranking <- ggplot(rank_wide, aes(MSE_rank, Cor_rank, color = Method)) +
  geom_vline(xintercept = mean(rank_wide$MSE_rank), linetype = "dotted", color = "grey55") +
  geom_hline(yintercept = mean(rank_wide$Cor_rank), linetype = "dotted", color = "grey55") +
  geom_point(size = 4) +
  geom_text_repel(aes(label = Method), size = 3.4, max.overlaps = 20, seed = 1) +
  scale_color_manual(values = PAL, drop = FALSE, guide = "none") +
  scale_x_reverse() + scale_y_reverse() +
  labs(title = "Accuracy vs ranking-quality tradeoff (avg rank across all cells)",
       subtitle = "Up & right = better on both. X: MSE magnitude (lower rank better). Y: correlation with true CATE (beta>0).",
       x = "Average MSE rank  (right = more accurate)",
       y = "Average correlation rank\n(up = better ranking)")

# (10) Pairwise win-rate heatmap among the CATE methods (from the Gao-2024
#      pairwise summary written by scripts/aggregate.R; skipped if absent)
if (file.exists("results/final_pairwise_relative_error.csv")) {
  pw <- read.csv("results/final_pairwise_relative_error.csv", check.names = FALSE)
  pw <- pw[pw$Method1 %in% FOCUS_CATE & pw$Method2 %in% FOCUS_CATE, ]
  # Most head-to-heads are statistical ties, so raw win% is small. Use the NET
  # win margin (row% - col%), diverging at 0: blue = row wins more often.
  win <- bind_rows(
    pw %>% transmute(Row = Method1, Col = Method2,
                     net = Pct_Sims_Method1_Wins - Pct_Sims_Method2_Wins),
    pw %>% transmute(Row = Method2, Col = Method1,
                     net = Pct_Sims_Method2_Wins - Pct_Sims_Method1_Wins)
  ) %>% group_by(Row, Col) %>% summarise(Net = mean(net), .groups = "drop")
  win$Row <- factor(win$Row, levels = FOCUS_CATE)
  win$Col <- factor(win$Col, levels = FOCUS_CATE)
  lim <- max(abs(win$Net))
  plots$Pairwise <- ggplot(win, aes(Col, fct_rev(Row), fill = Net)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%+.0f", Net)), size = 3) +
    scale_fill_gradient2(low = "#B2182B", mid = "#F7F7F7", high = "#2166AC",
                         midpoint = 0, limits = c(-lim, lim),
                         name = "Net win %\n(row - col)") +
    labs(title = "Pairwise head-to-head: net win margin (ROW minus COLUMN)",
         subtitle = "Avg over all scenarios & beta. Blue (+) = row beats column more often; most pairs tie in ~70-95% of sims.",
         x = "Column method", y = "Row method") +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
}

# (11) TMLE (Li 2026 VTE) — VTE estimate + reject rate vs beta (skipped if the
#      VTE method was not part of the run)
if (file.exists("results/final_vte_summary.csv")) {
  vte <- read.csv("results/final_vte_summary.csv", check.names = FALSE) %>%
    mutate(Scen = lab_scen(Scenario),
           lo = Mean_VTE - 1.96*SD_VTE/sqrt(N_Sims),
           hi = Mean_VTE + 1.96*SD_VTE/sqrt(N_Sims))
  plots$VTE_est <- ggplot(vte, aes(Beta, Mean_VTE, group = 1)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, fill = PAL["TMLE (Li 2026 VTE)"]) +
    geom_line(color = PAL["TMLE (Li 2026 VTE)"], linewidth = 0.9) +
    geom_point(color = PAL["TMLE (Li 2026 VTE)"], size = 1.9) +
    facet_wrap(~ Scen, scales = "free_y", ncol = 3) +
    labs(title = "TMLE (Li 2026 VTE): estimated variance of treatment effect",
         subtitle = "VTE grows with beta where heterogeneity is real. Ribbon = 95% MC CI.",
         x = expression(beta), y = "Mean VTE")
  plots$VTE_power <- ggplot(vte, aes(Beta, Reject_Rate_05, group = 1)) +
    geom_hline(yintercept = 5, linetype = "dashed", color = "grey40") +
    geom_line(color = PAL["TMLE (Li 2026 VTE)"], linewidth = 0.9) +
    geom_point(color = PAL["TMLE (Li 2026 VTE)"], size = 1.9) +
    facet_wrap(~ Scen, ncol = 3) +
    labs(title = "TMLE (Li 2026 VTE): rejection rate (heterogeneity test)",
         subtitle = "beta=0 = type-I error (target ~5%, dashed); beta>0 = power.",
         x = expression(beta), y = "Reject rate (%)")
}

# (12) Runtime
runt <- cell %>% group_by(Method) %>% summarise(Time = mean(Time), .groups = "drop")
plots$Runtime <- ggplot(runt, aes(reorder(Method, Time), Time, fill = Method)) +
  geom_col(alpha = 0.9) +
  geom_text(aes(label = sprintf("%.1fs", Time)), hjust = -0.15, size = 3) +
  scale_fill_manual(values = PAL, drop = FALSE) +
  coord_flip() + expand_limits(y = max(runt$Time)*1.15) +
  labs(title = "Mean runtime per fit", x = NULL, y = "Seconds (avg over all cells)") +
  theme(legend.position = "none")

# ---------------------------------------------------------------------
# 4. Emit: one multipage PDF + individual PNGs
# ---------------------------------------------------------------------
# Only plots that were actually built (guarded ones may be absent).
order <- intersect(
  c("MSE","Bias","Variance","Correlation","Somers","Detection",
    "DetectionVariants","MSE_dist","Heatmap","Ranking","Pairwise",
    "VTE_est","VTE_power","Runtime"),
  names(plots))
sizes <- list(default = c(13, 7.5), Pairwise = c(9, 7.5), Ranking = c(10, 7),
              Runtime = c(9, 5))

pdf(file.path(OUT, "Focus_Methods_Analysis.pdf"), width = 13, height = 7.5, onefile = TRUE)
for (nm in order) {
  print(plots[[nm]])
  sz <- if (!is.null(sizes[[nm]])) sizes[[nm]] else sizes$default
  ggsave(file.path(OUT, sprintf("%02d_%s.png", match(nm, order), nm)),
         plots[[nm]], width = sz[1], height = sz[2], dpi = 150)
}
invisible(dev.off())

# ---------------------------------------------------------------------
# 5. Self-contained report PDF: cover + findings + ranking table + plots
# ---------------------------------------------------------------------
methods_txt <- paste0("Methods analysed\n", paste0("  -  ", FOCUS, collapse = "\n"))
scen_txt <- paste(
  "Scenarios (CATE structure)",
  "  S1   Step in X11 (probit)",
  "  S2   Linear in X14",
  "  S3   Subgroup AND (X14 & X1)",
  "  S4   Subgroup OR (X14 | X4)", sep = "\n")
findings_txt <- paste(
  "Quantitative findings: read them directly from the ranking table on the next",
  "page and the numbered metric panels in this folder (focus_overall_ranking.csv,",
  "focus_summary_by_cell.csv). A hand-written narrative is intentionally omitted",
  "so the report cannot drift out of sync with the underlying data.",
  "",
  "Scenario set (Sun et al. 2022 benchtm): S1 step in X11 (probit), S2 linear in",
  sprintf("X14, S3 subgroup AND (X14 & X1), S4 subgroup OR (X14 | X4). beta in {%s}",
          paste(sort(unique(cell$Beta)), collapse = ",")),
  "(beta=0 = null / type-I error). True effect modifiers per scenario are derived",
  "from the DGP in R/scenario_meta.R::get_true_modifiers().",
  sep = "\n")
rank_tbl <- rank_wide %>% arrange(MSE_rank) %>%
  transmute(Method, `Avg MSE rank` = sprintf("%.2f", MSE_rank),
            `Avg correlation rank` = sprintf("%.2f", Cor_rank))

pdf(file.path(OUT, "Focus_Methods_Report.pdf"), width = 13, height = 7.5, onefile = TRUE)
# cover
grid.newpage()
grid.text("Focused HTE Meta-Learner Analysis", 0.5, 0.88, gp = gpar(fontsize = 26, fontface = "bold"))
grid.text(sprintf("Cross-scenario comparison of %d selected methods", length(FOCUS)),
          0.5, 0.80, gp = gpar(fontsize = 14, col = "grey30"))
grid.text(sprintf("Generated %s   |   RCT design   |   beta in {%s}   |   %s Monte-Carlo sims/cell",
                  Sys.Date(), paste(sort(unique(cell$Beta)), collapse = ", "),
                  paste(range(cell$n), collapse = "-")),
          0.5, 0.74, gp = gpar(fontsize = 11, col = "grey45"))
grid.text(methods_txt, 0.10, 0.60, just = c("left", "top"), gp = gpar(fontsize = 12))
grid.text(scen_txt,    0.58, 0.60, just = c("left", "top"), gp = gpar(fontsize = 12))
grid.text("beta = 0 is the null (no heterogeneity); larger beta = stronger HTE.",
          0.5, 0.08, gp = gpar(fontsize = 10, col = "grey45", fontface = "italic"))
# findings
grid.newpage()
grid.text("Key findings", 0.06, 0.95, just = c("left", "top"), gp = gpar(fontsize = 20, fontface = "bold"))
grid.text(findings_txt, 0.06, 0.86, just = c("left", "top"), gp = gpar(fontsize = 11, lineheight = 1.4))
# ranking table
grid.newpage()
grid.text("Overall average rank across all scenario x beta cells (1 = best)",
          0.5, 0.90, gp = gpar(fontsize = 16, fontface = "bold"))
grid.draw(tableGrob(rank_tbl, rows = NULL, theme = ttheme_minimal(base_size = 13)))
# all plots
for (nm in order) print(plots[[nm]])
invisible(dev.off())

cat("\nDone. Outputs in", OUT, ":\n")
cat(" - Focus_Methods_Report.pdf  (cover + findings + table + ", length(order), " plots )\n", sep = "")
cat(" - Focus_Methods_Analysis.pdf (", length(order), " plot pages )\n", sep = "")
cat(" - individual PNGs 01..", length(order), "\n", sep = "")
cat(" - focus_summary_by_cell.csv, focus_overall_ranking.csv\n")
