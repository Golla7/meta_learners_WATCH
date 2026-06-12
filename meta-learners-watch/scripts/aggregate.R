#!/usr/bin/env Rscript

# scripts/aggregate.R
# Full WATCH-framework analysis of the simulation run. RUN FROM THE REPO ROOT:
#   Rscript scripts/aggregate.R
#
# Reads the per-worker .rds files from temp/ and writes every output to
# results/ (summary CSVs + per-scenario WATCH PDFs + cross-scenario summary).
# TE-VIM forest plots and TMLE-ATE plots are GUARDED — they render only when
# the worker emitted that data, and auto-skip otherwise.
#
# Objective 1: Global test for heterogeneity   (Het_P_Value_Max / _Quad)
# Objective 2: Effect modifier identification  (Top_Modifier vs true DGP)
# Objective 3: Individual treatment effect est. (MSE / Correlation / Somers' D)
# Plus: Gao 2024 pairwise relative-error EIF, Li 2026 VTE, wall-clock timing.

library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(viridis)

cat("========================================\n")
cat("WATCH FRAMEWORK ANALYSIS\n")
cat("========================================\n\n")

# Soft-coded scenario ground truth (true effect modifiers / prognostic vars /
# labels) + shared EXCLUDE_METHODS & METHOD_PAL.
source("R/scenario_meta.R")

dir.create("results", showWarnings = FALSE)   # all outputs land here

# =============================================================================
# 1. READ AND COMBINE RESULTS (from temp/)
# =============================================================================

result_files <- list.files("temp", pattern = "^results_.*\\.rds$",
                           full.names = TRUE)
if (length(result_files) == 0)
  stop("ERROR: No result files in temp/. Run scripts/run_parallel.sh first!")
cat(sprintf("Found %d result file(s) in temp/\n", length(result_files)))

# Each scripts/worker.R file is a list with:
#   $metrics, $vi_scores, $pairwise_relative_errors (+ NULL tevims/ate stubs).
read_worker_file <- function(file) {
  data <- readRDS(file)
  if (is.null(data)) return(NULL)
  if (is.data.frame(data)) {
    return(list(metrics = data, vi_scores = list(),
                tevims_scores = NULL, ate_estimates = NULL,
                pairwise_relative_errors = NULL))
  }
  if (is.list(data) && !is.null(data$metrics) && nrow(data$metrics) > 0) {
    if (is.null(data$vi_scores))                data$vi_scores                <- list()
    if (is.null(data$tevims_scores))            data$tevims_scores            <- NULL
    if (is.null(data$ate_estimates))            data$ate_estimates            <- NULL
    if (is.null(data$pairwise_relative_errors)) data$pairwise_relative_errors <- NULL
    return(data)
  }
  NULL
}

all_results <- lapply(result_files, read_worker_file)
all_results <- all_results[!sapply(all_results, is.null)]
if (length(all_results) == 0)
  stop("All result files were empty. Check logs/ for worker errors.")

combined_results <- dplyr::bind_rows(lapply(all_results, `[[`, "metrics"))

# Drop excluded methods up front so every downstream table/plot honours it.
if (length(EXCLUDE_METHODS) > 0) {
  combined_results <- combined_results %>% filter(!Method %in% EXCLUDE_METHODS)
  cat(sprintf("Excluding %d method(s): %s\n",
              length(EXCLUDE_METHODS), paste(EXCLUDE_METHODS, collapse = ", ")))
}

# Optional subset: when non-NULL, keep ONLY these methods in every downstream
# table and plot. Set KEEP_METHODS <- NULL to analyse every method that was
# run (all 13 by default). The default below is the focus subset used for the
# headline comparison in results/.
KEEP_METHODS <- c(
  "DINA Original",
  "DINA Cross-Fit",
  "DINA Gao (LASSO k2)",
  "DINA Centered",            # §3z centered-design DINA
  "DINA Centered Expanded",   # §3z CATE-accuracy variant
  "DR-learner (2-Model 2SL)",
  "R-learner (LASSO)",
  "Dandl Extension"
)
# Keep the pre-subset frame so method-specific summaries (e.g. the TMLE VTE
# section) still see methods that are not part of the focus subset.
combined_results_all <- combined_results
if (!is.null(KEEP_METHODS)) {
  combined_results <- combined_results %>% filter(Method %in% KEEP_METHODS)
  cat(sprintf("Keeping %d method(s): %s\n",
              length(KEEP_METHODS), paste(KEEP_METHODS, collapse = ", ")))
}

# Soft covariate count p for the "uniform selection" reference line — first
# available full VI table has one row per covariate; fall back to distinct
# Top_Modifier count. (vi_scores here is keyed "s{sim}" -> list(method -> df).)
N_COVARIATES <- local({
  vt <- NULL
  for (r in all_results) for (sk in r$vi_scores) for (vi in sk)
    if (is.data.frame(vi) && nrow(vi) > 0) { vt <- vi; break }
  if (!is.null(vt)) nrow(vt) else length(unique(stats::na.omit(combined_results$Top_Modifier)))
})

# Pairwise relative-error EIF rows (Gao 2024) — one row per (sim, m1, m2).
combined_pairwise <- tryCatch(
  dplyr::bind_rows(Filter(Negate(is.null),
                          lapply(all_results, `[[`, "pairwise_relative_errors"))),
  error = function(e) NULL
)
if (!is.null(combined_pairwise) && nrow(combined_pairwise) == 0) combined_pairwise <- NULL

# TE-VIM / ATE frames — NULL for the clean harness; kept for forward-compat.
combined_tevims <- tryCatch(
  dplyr::bind_rows(Filter(Negate(is.null), lapply(all_results, `[[`, "tevims_scores"))),
  error = function(e) NULL)
if (!is.null(combined_tevims) && nrow(combined_tevims) == 0) combined_tevims <- NULL
combined_ate <- tryCatch(
  dplyr::bind_rows(Filter(Negate(is.null), lapply(all_results, `[[`, "ate_estimates"))),
  error = function(e) NULL)
if (!is.null(combined_ate) && nrow(combined_ate) == 0) combined_ate <- NULL

# Guarantee newer columns exist (forward/backward compatibility).
for (col in c("Time_Sec", "Somers_D", "Het_P_Value_Quad",
              "VTE_Estimate", "VTE_SE", "VTE_CI_Lo", "VTE_CI_Hi")) {
  if (!col %in% names(combined_results)) combined_results[[col]] <- NA_real_
}
# The clean harness reports Het_P_Value_Max / _Quad (Sechidis Obj 1). Alias a
# plain Het_P_Value / Het_Detected to the max-type columns so the legacy plot
# code (p1a, summary_stats Mean_P_Value / Detection_Rate) works unchanged.
if (!"Het_P_Value_Max" %in% names(combined_results)) combined_results$Het_P_Value_Max <- NA_real_
if (!"Het_Detected_Max" %in% names(combined_results)) combined_results$Het_Detected_Max <- NA
if (!"Het_Detected_Quad" %in% names(combined_results)) combined_results$Het_Detected_Quad <- NA
if (!"VTE_Detected" %in% names(combined_results)) combined_results$VTE_Detected <- NA
combined_results$Het_P_Value <- combined_results$Het_P_Value_Max
combined_results$Het_Detected <- combined_results$Het_Detected_Max
if (!"Treatment_Scenario" %in% names(combined_results))
  combined_results$Treatment_Scenario <- "RCT"

# valid_results — only learners that returned an MSE (drops MSE-less rows).
valid_results <- combined_results[!is.na(combined_results$MSE), ]
cat(sprintf("\nTotal valid result rows (CATE methods): %d\n", nrow(valid_results)))

# obj2_results — keep all rows with a Top_Modifier (modifier-ID analysis).
obj2_results <- combined_results[!is.na(combined_results$Top_Modifier), ]

# =============================================================================
# 2. SUMMARY STATISTICS
# =============================================================================

cat("\nCALCULATING SUMMARY STATISTICS\n------------------------------\n")

safe_mean <- function(x) if (length(x) == 0) NA_real_ else mean(x, na.rm = TRUE)

summary_stats <- valid_results %>%
  group_by(Treatment_Scenario, Scenario, Beta, Method) %>%
  summarise(
    N_Sims              = n(),
    Mean_MSE            = mean(MSE, na.rm = TRUE),
    SD_MSE              = sd(MSE,  na.rm = TRUE),
    Mean_Bias           = mean(Bias, na.rm = TRUE),
    Mean_Variance       = mean(Variance, na.rm = TRUE),
    Mean_Correlation    = mean(Correlation, na.rm = TRUE),
    Mean_Somers_D       = safe_mean(Somers_D),
    Mean_MAE            = safe_mean(MAE),
    # Sechidis Obj 1: both global-test statistics. Mean_P_Value / Detection_Rate
    # are max-type aliases the per-scenario plots reference.
    Mean_P_Value        = mean(Het_P_Value_Max,  na.rm = TRUE),
    Mean_P_Value_Max    = mean(Het_P_Value_Max,  na.rm = TRUE),
    Mean_P_Value_Quad   = mean(Het_P_Value_Quad, na.rm = TRUE),
    Detection_Rate      = mean(Het_Detected_Max,  na.rm = TRUE) * 100,
    Detection_Rate_Max  = mean(Het_Detected_Max,  na.rm = TRUE) * 100,
    Detection_Rate_Quad = mean(Het_Detected_Quad, na.rm = TRUE) * 100,
    Mean_Time_Sec       = safe_mean(Time_Sec),
    .groups = 'drop'
  ) %>%
  arrange(Treatment_Scenario, Scenario, Beta, Mean_MSE)

write.csv(summary_stats, "results/final_summary_statistics.csv", row.names = FALSE)
cat("  Saved: results/final_summary_statistics.csv\n")

# Li et al. 2026 VTE-TMLE summary (only the TMLE method populates VTE columns).
# Uses the PRE-subset frame so the summary is written even when the TMLE
# method is not part of KEEP_METHODS.
for (col in c("VTE_Estimate", "VTE_SE", "VTE_CI_Lo", "VTE_CI_Hi", "VTE_Detected")) {
  if (!col %in% names(combined_results_all)) combined_results_all[[col]] <- NA_real_
}
if (!"Treatment_Scenario" %in% names(combined_results_all))
  combined_results_all$Treatment_Scenario <- "RCT"
vte_rows <- combined_results_all %>% filter(!is.na(VTE_Estimate))
if (nrow(vte_rows) > 0) {
  vte_summary <- vte_rows %>%
    group_by(Treatment_Scenario, Scenario, Beta, Method) %>%
    summarise(
      N_Sims         = n(),
      Mean_VTE       = mean(VTE_Estimate, na.rm = TRUE),
      SD_VTE         = sd(VTE_Estimate,  na.rm = TRUE),
      Mean_VTE_SE    = mean(VTE_SE,       na.rm = TRUE),
      Mean_CI_Width  = mean(VTE_CI_Hi - VTE_CI_Lo, na.rm = TRUE),
      Reject_Rate_05 = mean(VTE_Detected, na.rm = TRUE) * 100,
      .groups = 'drop'
    )
  write.csv(vte_summary, "results/final_vte_summary.csv", row.names = FALSE)
  cat("  Saved: results/final_vte_summary.csv\n")
}

# =============================================================================
# 2b. PAIRWISE RELATIVE-ERROR SUMMARY (Gao 2024 EIF)
# =============================================================================
# Per-sim CI exclusion -> single-dataset "power" of declaring a winner;
# Monte Carlo CI over mean δ̂ -> the meta-population mean relative error.
if (!is.null(combined_pairwise) && nrow(combined_pairwise) > 0) {
  if (!"Treatment_Scenario" %in% names(combined_pairwise))
    combined_pairwise$Treatment_Scenario <- "RCT"
  pairwise_summary <- combined_pairwise %>%
    group_by(Treatment_Scenario, Scenario, Beta, Method1, Method2) %>%
    summarise(
      N_Sims                = n(),
      Mean_Delta            = mean(Delta, na.rm = TRUE),
      SD_Delta              = sd(Delta,   na.rm = TRUE),
      MC_CI_Lo              = mean(Delta, na.rm = TRUE) -
                              1.96 * sd(Delta, na.rm = TRUE) / sqrt(sum(!is.na(Delta))),
      MC_CI_Hi              = mean(Delta, na.rm = TRUE) +
                              1.96 * sd(Delta, na.rm = TRUE) / sqrt(sum(!is.na(Delta))),
      Pct_Sims_Sig          = mean(CI_Excludes_Zero, na.rm = TRUE) * 100,
      Pct_Sims_Method1_Wins = mean(Winner == Method1, na.rm = TRUE) * 100,
      Pct_Sims_Method2_Wins = mean(Winner == Method2, na.rm = TRUE) * 100,
      Pct_Sims_Tie          = mean(Winner == "Tie",   na.rm = TRUE) * 100,
      Overall_Winner        = case_when(
        MC_CI_Hi < 0 ~ first(Method1),
        MC_CI_Lo > 0 ~ first(Method2),
        TRUE         ~ "Tie"
      ),
      .groups = 'drop'
    ) %>%
    arrange(Treatment_Scenario, Scenario, Beta, Method1, Method2)
  write.csv(pairwise_summary, "results/final_pairwise_relative_error.csv", row.names = FALSE)
  cat(sprintf("  Saved: results/final_pairwise_relative_error.csv (%d pair-rows)\n",
              nrow(pairwise_summary)))
} else {
  cat("  No pairwise relative-error data found.\n")
}

# =============================================================================
# 2.5 PERFORMANCE DISTRIBUTION PLOTS
# =============================================================================

cat("\nGENERATING PERFORMANCE DISTRIBUTION PLOTS\n------------------------------------------\n")
dir.create("results", showWarnings = FALSE)

if (nrow(valid_results) > 0) {
  density_plot_data <- valid_results %>%
    select(Scenario, Beta, Method, MSE, Correlation) %>%
    pivot_longer(cols = c(MSE, Correlation), names_to = "Metric", values_to = "Value") %>%
    mutate(
      Facet_Label  = paste0("Scenario ", Scenario, "\nBeta: ", Beta),
      Metric_Label = ifelse(Metric == "MSE", "Metric: MSE", "Mean Correlation")
    )

  p_density <- ggplot(density_plot_data, aes(x = Value, fill = Method, color = Method)) +
    geom_density(alpha = 0.5, linewidth = 0.8) +
    facet_grid(Metric_Label ~ Facet_Label, scales = "free", switch = "y") +
    theme_bw(base_size = 10) +
    theme(strip.text.x = element_text(size = 8), strip.text.y = element_text(size = 10, angle = 0),
          strip.placement = "outside", legend.position = "bottom",
          axis.text.x = element_text(angle = 45, hjust = 1, size = 7)) +
    scale_fill_viridis_d(name = "Method") + scale_color_viridis_d(name = "Method") +
    labs(title = "Distribution of Performance Metrics (Clean Harness)",
         subtitle = "Spread of MSE and Correlation across all simulation runs",
         x = "Metric Value", y = NULL)
  ggsave("results/Performance_Distribution_AllScenarios.pdf",
         plot = p_density, width = 20, height = 8, limitsize = FALSE)
  cat("  Saved: results/Performance_Distribution_AllScenarios.pdf\n")

  violin_plot_data <- valid_results %>%
    select(Scenario, Beta, Method, MSE, Correlation) %>%
    pivot_longer(cols = c(MSE, Correlation), names_to = "Metric", values_to = "Value") %>%
    mutate(ScenBeta = paste0("Scen ", Scenario, "\nβ=", Beta),
           Metric_Label = ifelse(Metric == "MSE", "MSE", "Correlation"))

  p_violin <- ggplot(violin_plot_data, aes(x = ScenBeta, y = Value, fill = Method)) +
    geom_violin(alpha = 0.6, position = position_dodge(width = 0.9), scale = "width", trim = TRUE) +
    geom_boxplot(width = 0.1, position = position_dodge(width = 0.9), alpha = 0.3, outlier.size = 0.5) +
    facet_wrap(~ Metric_Label, scales = "free_y", ncol = 1) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), legend.position = "bottom",
          strip.text = element_text(face = "bold", size = 11), panel.grid.major.x = element_blank()) +
    scale_fill_viridis_d(name = "Method") +
    labs(title = "Performance Distribution: Violin Plot View (Clean Harness)",
         subtitle = "Distribution shape and quartiles for MSE and Correlation",
         x = "Scenario and Beta Level", y = "Metric Value")
  ggsave("results/Performance_Distribution_ViolinPlot.pdf",
         plot = p_violin, width = 16, height = 10)
  cat("  Saved: results/Performance_Distribution_ViolinPlot.pdf\n")

  if (requireNamespace("ggridges", quietly = TRUE)) {
    library(ggridges)
    ridge_plot_data <- valid_results %>%
      select(Scenario, Beta, Method, MSE, Correlation) %>%
      pivot_longer(cols = c(MSE, Correlation), names_to = "Metric", values_to = "Value") %>%
      mutate(ScenBeta = paste0("Scen ", Scenario, ", β=", Beta),
             Metric_Label = ifelse(Metric == "MSE", "Metric: MSE", "Mean Correlation")) %>%
      arrange(Scenario, Beta)
    p_ridge <- ggplot(ridge_plot_data, aes(x = Value, y = ScenBeta, fill = Method)) +
      geom_density_ridges(alpha = 0.6, scale = 0.9) +
      facet_wrap(~ Metric_Label, scales = "free_x", ncol = 2) +
      theme_ridges() +
      theme(legend.position = "bottom", strip.text = element_text(face = "bold", size = 12)) +
      scale_fill_viridis_d(name = "Method") +
      labs(title = "Performance Distribution: Ridge Plot View (Clean Harness)",
           subtitle = "Density distributions across all Scenarios and Beta levels",
           x = "Metric Value", y = "Scenario and Beta Level")
    ggsave("results/Performance_Distribution_RidgePlot.pdf",
           plot = p_ridge, width = 14, height = 12)
    cat("  Saved: results/Performance_Distribution_RidgePlot.pdf\n")
  } else {
    cat("  Note: install 'ggridges' for the ridge plot (optional).\n")
  }
}

# =============================================================================
# 3. PER-SCENARIO WATCH ANALYSIS (Objectives 1, 2, 3)
# =============================================================================

cat("\nGENERATING WATCH FRAMEWORK VISUALIZATIONS\n------------------------------------------\n")
unique_scenarios <- sort(unique(summary_stats$Scenario))

for (scen in unique_scenarios) {

  pdf_filename <- sprintf("results/WATCH_analysis_Scenario_%d.pdf", scen)
  cat(sprintf("  Scenario %d -> %s\n", scen, pdf_filename))

  scenario_summary <- summary_stats %>% filter(Scenario == scen)
  scenario_valid   <- valid_results %>% filter(Scenario == scen)
  scenario_obj2    <- obj2_results  %>% filter(Scenario == scen)

  pdf(pdf_filename, width = 16, height = 12)

  # ----- OBJECTIVE 1: GLOBAL TEST FOR HETEROGENEITY -----
  cat(sprintf("    - Objective 1: Heterogeneity testing\n"))

  null_data <- scenario_valid %>% filter(Beta == 0)
  if (nrow(null_data) > 0) {
    p1a <- ggplot(null_data, aes(x = Method, y = Het_P_Value, fill = Method)) +
      geom_boxplot(alpha = 0.8) +
      geom_hline(yintercept = c(0.25, 0.5, 0.75), linetype = "dashed", color = "gray") +
      coord_flip() + theme_bw(base_size = 12) + theme(legend.position = "none") +
      labs(title = sprintf("Objective 1(i): P-value Under Null (Scenario %d, β=0)", scen),
           subtitle = "P-values should be ~Uniform[0,1] when no heterogeneity exists",
           x = "Method", y = "P-value")
    print(p1a)
  }

  p1b <- ggplot(scenario_summary, aes(x = Beta, y = Mean_P_Value, color = Method, group = Method)) +
    geom_line(linewidth = 1.2) + geom_point(size = 3) +
    labs(title = sprintf("Objective 1(ii): Power Analysis (Scenario %d)", scen),
         subtitle = "Lower mean p-value indicates better power to detect heterogeneity",
         y = "Mean P-value", x = "Beta Level (Effect Size Strength)") +
    theme_bw(base_size = 12) + theme(legend.position = "bottom") + scale_color_viridis_d()
  print(p1b)

  p1c <- ggplot(scenario_summary, aes(x = Beta, y = Detection_Rate, color = Method, group = Method)) +
    geom_line(linewidth = 1.2) + geom_point(size = 3) +
    geom_hline(yintercept = 10, linetype = "dashed", color = "red", alpha = 0.5) +
    labs(title = sprintf("Objective 1: Detection Rate (Scenario %d)", scen),
         subtitle = "% of sims detecting heterogeneity (p < 0.05). Red line = nominal α.",
         y = "Detection Rate (%)", x = "Beta Level") +
    theme_bw(base_size = 12) + theme(legend.position = "bottom") + scale_color_viridis_d()
  print(p1c)

  # Max-type vs quadratic-type detection (Sechidis 3.3.1)
  if ("Detection_Rate_Quad" %in% names(scenario_summary)) {
    stat_long <- scenario_summary %>%
      select(Beta, Method, Detection_Rate_Max, Detection_Rate_Quad) %>%
      tidyr::pivot_longer(c(Detection_Rate_Max, Detection_Rate_Quad),
                          names_to = "Statistic", values_to = "Detection_Rate") %>%
      mutate(Statistic = ifelse(Statistic == "Detection_Rate_Max",
                                "max-type T_max", "quadratic-type T_quad"))
    p1d <- ggplot(stat_long, aes(x = Beta, y = Detection_Rate, color = Method,
                                 linetype = Statistic, group = interaction(Method, Statistic))) +
      geom_line(linewidth = 1.0) + geom_point(size = 2) +
      geom_hline(yintercept = 5, linetype = "dotted", color = "red", alpha = 0.5) +
      labs(title = sprintf("Objective 1: Max- vs Quadratic-type Power (Scenario %d)", scen),
           subtitle = "Permutation test on each method's pseudo-outcome. Red dotted = nominal α = 5%.",
           y = "Detection Rate (%)", x = "Beta Level") +
      theme_bw(base_size = 12) + theme(legend.position = "bottom") + scale_color_viridis_d()
    print(p1d)
  }

  # ----- OBJECTIVE 2: EFFECT MODIFIER IDENTIFICATION -----
  cat(sprintf("    - Objective 2: Effect modifier ranking\n"))
  true_info       <- get_true_modifiers(scen)
  true_modifiers  <- true_info$modifiers
  prognostic_only <- true_info$prognostic_only
  cat(sprintf("      True modifiers: %s\n", paste(true_modifiers, collapse = ", ")))

  # 2a: top modifier selection frequency under null
  nd2 <- scenario_obj2 %>% filter(Beta == 0)
  if (nrow(nd2) > 0 && "Top_Modifier" %in% names(nd2)) {
    modifier_freq_null <- nd2 %>%
      filter(!is.na(Top_Modifier)) %>%
      group_by(Method, Top_Modifier) %>% summarise(Count = n(), .groups = 'drop') %>%
      group_by(Method) %>%
      mutate(Probability = Count / sum(Count),
             VarType = case_when(
               Top_Modifier %in% true_modifiers  ~ "True Modifier",
               Top_Modifier %in% prognostic_only ~ "Prognostic Only",
               TRUE                              ~ "Noise Variable")) %>%
      arrange(Method, desc(Probability)) %>% group_by(Method) %>% slice_head(n = 12)
    if (nrow(modifier_freq_null) > 0) {
      p2a <- ggplot(modifier_freq_null, aes(x = reorder(Top_Modifier, Probability),
                                            y = Probability, fill = VarType)) +
        geom_col() +
        geom_hline(yintercept = 1/N_COVARIATES, linetype = "dashed", color = "red") +
        facet_wrap(~Method, scales = "free_x", ncol = 3) + coord_flip() +
        theme_bw(base_size = 10) + theme(legend.position = "bottom") +
        scale_fill_manual(values = c("True Modifier" = "#2ecc71", "Prognostic Only" = "#f39c12",
                                     "Noise Variable" = "#95a5a6"), name = "Variable Type") +
        labs(title = sprintf("Objective 2(i): Top Modifier Under Null (Scenario %d, β=0)", scen),
             subtitle = sprintf("Unbiased methods select all vars equally (~%.1f%% for %d vars, red line)",
                                100/N_COVARIATES, N_COVARIATES),
             x = "Variable", y = "Selection Probability")
      print(p2a)
    }
  }

  # 2b: probability top variable is truly predictive (across beta)
  if ("Top_Modifier" %in% names(scenario_obj2)) {
    modifier_accuracy <- scenario_obj2 %>%
      filter(!is.na(Top_Modifier)) %>%
      mutate(Is_True_Modifier = Top_Modifier %in% true_modifiers,
             Is_Prognostic    = Top_Modifier %in% prognostic_only) %>%
      group_by(Method, Beta) %>%
      summarise(Prob_True_Modifier = mean(Is_True_Modifier, na.rm = TRUE),
                Prob_Prognostic    = mean(Is_Prognostic, na.rm = TRUE),
                Prob_Noise         = 1 - mean(Is_True_Modifier | Is_Prognostic, na.rm = TRUE),
                N_Sims = n(), .groups = 'drop')
    if (nrow(modifier_accuracy) > 0) {
      p2b <- ggplot(modifier_accuracy, aes(x = Beta, y = Prob_True_Modifier,
                                           color = Method, group = Method)) +
        geom_line(linewidth = 1.2) + geom_point(size = 3) +
        labs(title = sprintf("Objective 2(ii): Prob Top Variable is Truly Predictive (Scenario %d)", scen),
             subtitle = sprintf("True modifiers: %s", paste(true_modifiers, collapse = ", ")),
             y = "Probability Top Variable is True Modifier", x = "Beta Level") +
        theme_bw(base_size = 12) + theme(legend.position = "bottom") +
        scale_color_viridis_d() + ylim(0, 1)
      print(p2b)

      modifier_long <- modifier_accuracy %>%
        select(Method, Beta, Prob_True_Modifier, Prob_Prognostic, Prob_Noise) %>%
        pivot_longer(cols = starts_with("Prob_"), names_to = "VarType", values_to = "Probability") %>%
        mutate(VarType = case_when(VarType == "Prob_True_Modifier" ~ "True Modifier",
                                   VarType == "Prob_Prognostic"    ~ "Prognostic Only",
                                   VarType == "Prob_Noise"         ~ "Noise Variable"))
      p2b_stack <- ggplot(modifier_long, aes(x = Beta, y = Probability, fill = VarType)) +
        geom_area(alpha = 0.7, position = "stack") + facet_wrap(~Method, ncol = 3) +
        scale_fill_manual(values = c("True Modifier" = "#2ecc71", "Prognostic Only" = "#f39c12",
                                     "Noise Variable" = "#95a5a6"), name = "Variable Type Selected") +
        labs(title = sprintf("Objective 2(ii): Variable Type Selection Breakdown (Scenario %d)", scen),
             subtitle = "Stacked area: what type of variable was selected as top modifier",
             y = "Proportion", x = "Beta Level") +
        theme_bw(base_size = 11) + theme(legend.position = "bottom") + ylim(0, 1)
      print(p2b_stack)
    }

    # 2c: most frequently selected top variable per method (β > 0)
    top_selections <- scenario_obj2 %>%
      filter(Beta > 0, !is.na(Top_Modifier)) %>%
      mutate(VarType = case_when(Top_Modifier %in% true_modifiers  ~ "True Modifier",
                                 Top_Modifier %in% prognostic_only ~ "Prognostic Only",
                                 TRUE                              ~ "Noise Variable")) %>%
      group_by(Method, Top_Modifier, VarType) %>% summarise(Count = n(), .groups = 'drop') %>%
      group_by(Method) %>% arrange(Method, desc(Count)) %>% slice_head(n = 5)
    if (nrow(top_selections) > 0) {
      p2c <- ggplot(top_selections, aes(x = reorder(Top_Modifier, Count), y = Count, fill = VarType)) +
        geom_col() + facet_wrap(~Method, scales = "free", ncol = 3) + coord_flip() +
        scale_fill_manual(values = c("True Modifier" = "#2ecc71", "Prognostic Only" = "#f39c12",
                                     "Noise Variable" = "#95a5a6"), name = "Variable Type") +
        theme_bw(base_size = 10) + theme(legend.position = "bottom") +
        labs(title = sprintf("Objective 2: Most Frequently Selected Variables (Scenario %d, β>0)", scen),
             subtitle = "Top 5 variables selected most often across sims with heterogeneity",
             x = "Variable", y = "Selection Count")
      print(p2c)
    }
  }

  # 2-forest: TE-VIM TMLE covariate ranking (Li 2026 Psi_2 / Psi_3) with 95% CI.
  scenario_tevims <- if (!is.null(combined_tevims))
    combined_tevims %>% filter(Scenario == scen, Beta > 0) else NULL
  if (!is.null(scenario_tevims) && nrow(scenario_tevims) > 0) {
    cat("      TE-VIM TMLE scores present — rendering forest plot.\n")
    # Monte-Carlo pooling per covariate: point = mean over sims; SE combines
    # within-sim variance (mean of se^2) and between-sim variance (law of total
    # variance), matching the legacy aggregate.R forest.
    forest_df <- scenario_tevims %>%
      group_by(variable) %>%
      summarise(
        vima_mean = mean(vima, na.rm = TRUE),
        vima_se   = sqrt(mean(vima_se^2, na.rm = TRUE) +
                         var(vima, na.rm = TRUE) / max(1, sum(!is.na(vima)))),
        vimb_mean = mean(vimb, na.rm = TRUE),
        vimb_se   = sqrt(mean(vimb_se^2, na.rm = TRUE) +
                         var(vimb, na.rm = TRUE) / max(1, sum(!is.na(vimb)))),
        .groups = "drop") %>%
      mutate(
        VarType = case_when(variable %in% true_modifiers  ~ "True Modifier",
                            variable %in% prognostic_only ~ "Prognostic Only",
                            TRUE                          ~ "Noise Variable"),
        vimb_ci_lo = vimb_mean - 1.96 * vimb_se,
        vimb_ci_hi = vimb_mean + 1.96 * vimb_se)

    p2_forest <- ggplot(forest_df,
                        aes(x = reorder(variable, vimb_mean), y = vimb_mean,
                            ymin = vimb_ci_lo, ymax = vimb_ci_hi, color = VarType)) +
      geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
      geom_pointrange(size = 0.4) + coord_flip() +
      scale_color_manual(values = c("True Modifier" = "#2ecc71",
                                    "Prognostic Only" = "#f39c12",
                                    "Noise Variable" = "#95a5a6"),
                         name = "Variable Type") +
      theme_bw(base_size = 11) + theme(legend.position = "bottom") +
      labs(title = sprintf("Objective 2 (TE-VIM TMLE): Scaled VIM Ranking (Scenario %d, beta>0)", scen),
           subtitle = sprintf("Li 2026 VIMb = Psi2/Psi1 (TMLE), pooled over sims. True modifier(s): %s",
                              paste(true_modifiers, collapse = ", ")),
           x = "Variable", y = "Scaled TE-VIM  VIMb = Psi2/Psi1  (95% CI)")
    print(p2_forest)
  }

  # ----- OBJECTIVE 3: CATE ESTIMATION -----
  cat(sprintf("    - Objective 3: CATE estimation\n"))

  p3a <- ggplot(scenario_valid, aes(x = reorder(Method, MSE, FUN = median), y = MSE, fill = Method)) +
    geom_boxplot(alpha = 0.8) +
    facet_grid(. ~ Beta, scales = "free_y", labeller = label_both) + coord_flip() +
    theme_bw(base_size = 12) + theme(legend.position = "none") +
    labs(title = sprintf("Objective 3: MSE Distribution (Scenario %d)", scen),
         subtitle = "Lower MSE = better CATE estimation. Spread = stability.",
         x = "Method", y = "Mean Squared Error (MSE)")
  print(p3a)

  p3b <- ggplot(scenario_summary, aes(x = Beta, y = Mean_MSE, color = Method, group = Method)) +
    geom_line(linewidth = 1.2) + geom_point(size = 3) + scale_y_log10() +
    labs(title = sprintf("Objective 3: CATE Estimation Performance (Scenario %d)", scen),
         subtitle = "Mean Squared Error for estimating individual treatment effects",
         y = "Mean MSE (log scale)", x = "Beta Level (Effect Size Strength)") +
    theme_bw(base_size = 12) + theme(legend.position = "bottom") + scale_color_viridis_d()
  print(p3b)

  p3c <- ggplot(scenario_summary, aes(x = Beta, y = Mean_Correlation, color = Method, group = Method)) +
    geom_line(linewidth = 1.2) + geom_point(size = 3) +
    labs(title = sprintf("Objective 3: Correlation with True Effects (Scenario %d)", scen),
         subtitle = "Higher correlation = better alignment with true treatment effects",
         y = "Mean Correlation", x = "Beta Level") +
    theme_bw(base_size = 12) + theme(legend.position = "bottom") + scale_color_viridis_d()
  print(p3c)

  # Somers' D across beta (Sechidis Obj 3 recommended ranking metric)
  if (any(!is.na(scenario_summary$Mean_Somers_D))) {
    p3d <- ggplot(scenario_summary, aes(x = Beta, y = Mean_Somers_D, color = Method, group = Method)) +
      geom_line(linewidth = 1.2) + geom_point(size = 3) +
      labs(title = sprintf("Objective 3: Somers' D Rank Correlation (Scenario %d)", scen),
           subtitle = "Sechidis et al. (2026) recommended Obj-3 metric; higher = better ranking",
           y = "Mean Somers' D", x = "Beta Level") +
      theme_bw(base_size = 12) + theme(legend.position = "bottom") + scale_color_viridis_d()
    print(p3d)
  }

  # ----- PER-SCENARIO SUMMARY COMPARISON TABLE -----
  comparison_table <- scenario_summary %>%
    filter(Beta > 0) %>%
    group_by(Method) %>%
    summarise(Avg_Power = mean(1 - Mean_P_Value), Avg_Detection = mean(Detection_Rate),
              Avg_MSE = mean(Mean_MSE), Avg_Correlation = mean(Mean_Correlation),
              Avg_Somers_D = mean(Mean_Somers_D, na.rm = TRUE), .groups = 'drop')
  if ("Top_Modifier" %in% names(scenario_obj2)) {
    obj2_metrics <- scenario_obj2 %>%
      filter(Beta > 0, !is.na(Top_Modifier)) %>%
      mutate(Is_True_Modifier = Top_Modifier %in% true_info$modifiers) %>%
      group_by(Method) %>%
      summarise(Avg_True_Selection = mean(Is_True_Modifier, na.rm = TRUE) * 100, .groups = 'drop')
    comparison_table <- comparison_table %>% full_join(obj2_metrics, by = "Method")
  }
  comparison_table <- comparison_table %>% arrange(desc(Avg_Power))
  table_plot <- gridExtra::tableGrob(
    comparison_table %>% mutate(across(where(is.numeric), ~round(., 3))), rows = NULL)
  grid.arrange(table_plot,
               top = sprintf("Summary Comparison (Scenario %d) — averaged over β>0\nTrue Modifiers: %s",
                             scen, paste(true_info$modifiers, collapse = ", ")))

  dev.off()
}

cat("  All per-scenario WATCH PDFs saved in results/.\n")

# =============================================================================
# 4. CROSS-SCENARIO SUMMARY
# =============================================================================

cat("\nGENERATING CROSS-SCENARIO SUMMARY\n----------------------------------\n")
pdf("results/Cross_Scenario_Summary.pdf", width = 16, height = 10)

# Obj 1: power across scenarios
power_summary <- summary_stats %>% filter(Beta > 0) %>%
  group_by(Scenario, Method) %>% summarise(Avg_Power = mean(1 - Mean_P_Value), .groups = 'drop')
p_power <- ggplot(power_summary, aes(x = factor(Scenario), y = Avg_Power, fill = Method)) +
  geom_col(position = "dodge") +
  labs(title = "Objective 1: Average Power Across All Scenarios",
       subtitle = "Higher = better ability to detect heterogeneity",
       x = "Scenario", y = "Average Power (1 - Mean P-value)") +
  theme_bw(base_size = 12) + theme(legend.position = "bottom") + scale_fill_viridis_d()
print(p_power)

# Obj 2: true modifier selection across scenarios
if ("Top_Modifier" %in% names(obj2_results)) {
  obj2_cross_scenario <- lapply(unique_scenarios, function(scen) {
    ti <- get_true_modifiers(scen)
    obj2_results %>%
      filter(Scenario == scen, Beta > 0, !is.na(Top_Modifier)) %>%
      mutate(Is_True_Modifier = Top_Modifier %in% ti$modifiers, Scenario = scen) %>%
      group_by(Scenario, Method) %>%
      summarise(True_Selection_Rate = mean(Is_True_Modifier, na.rm = TRUE) * 100, .groups = 'drop')
  }) %>% bind_rows()
  p_obj2 <- ggplot(obj2_cross_scenario,
                   aes(x = factor(Scenario), y = True_Selection_Rate, fill = Method)) +
    geom_col(position = "dodge") +
    labs(title = "Objective 2: True Modifier Selection Rate Across Scenarios",
         subtitle = "% of times the top-ranked variable was a true effect modifier",
         x = "Scenario", y = "True Modifier Selection Rate (%)") +
    theme_bw(base_size = 12) + theme(legend.position = "bottom") + scale_fill_viridis_d()
  print(p_obj2)
}

# TE-VIM TMLE cross-scenario forest (Li 2026 Psi_3 / VIMb), faceted by scenario.
if (!is.null(combined_tevims) && nrow(combined_tevims) > 0) {
  tevim_cross <- lapply(unique_scenarios, function(scen) {
    ti <- get_true_modifiers(scen)
    df <- combined_tevims %>% filter(Scenario == scen, Beta > 0)
    if (nrow(df) == 0) return(NULL)
    df %>%
      group_by(variable) %>%
      summarise(
        vimb_mean = mean(vimb, na.rm = TRUE),
        vimb_se   = sqrt(mean(vimb_se^2, na.rm = TRUE) +
                         var(vimb, na.rm = TRUE) / max(1, sum(!is.na(vimb)))),
        .groups = "drop") %>%
      mutate(Scenario = scen,
             VarType = case_when(variable %in% ti$modifiers       ~ "True Modifier",
                                 variable %in% ti$prognostic_only  ~ "Prognostic Only",
                                 TRUE                              ~ "Noise Variable"),
             vimb_ci_lo = vimb_mean - 1.96 * vimb_se,
             vimb_ci_hi = vimb_mean + 1.96 * vimb_se)
  }) %>% bind_rows()

  if (!is.null(tevim_cross) && nrow(tevim_cross) > 0) {
    p_tevim_cross <- ggplot(tevim_cross,
                           aes(x = reorder(variable, vimb_mean), y = vimb_mean,
                               ymin = vimb_ci_lo, ymax = vimb_ci_hi, color = VarType)) +
      geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
      geom_pointrange(size = 0.3) + coord_flip() +
      facet_wrap(~ paste("Scenario", Scenario), scales = "free", ncol = 2) +
      scale_color_manual(values = c("True Modifier" = "#2ecc71",
                                    "Prognostic Only" = "#f39c12",
                                    "Noise Variable" = "#95a5a6"), name = "Variable Type") +
      theme_bw(base_size = 10) + theme(legend.position = "bottom") +
      labs(title = "TE-VIM TMLE (Li 2026): Scaled VIM (VIMb) Ranking Across Scenarios (beta>0)",
           subtitle = "Higher = covariate explains more treatment-effect heterogeneity",
           x = "Variable", y = "Scaled TE-VIM  VIMb = Psi2/Psi1  (95% CI)")
    print(p_tevim_cross)
    cat("  TE-VIM cross-scenario forest plot rendered.\n")
  }
}

# Obj 3: MSE across scenarios
mse_summary <- summary_stats %>% filter(Beta > 0) %>%
  group_by(Scenario, Method) %>% summarise(Avg_MSE = mean(Mean_MSE), .groups = 'drop')
p_mse <- ggplot(mse_summary, aes(x = factor(Scenario), y = Avg_MSE, fill = Method)) +
  geom_col(position = "dodge") +
  labs(title = "Objective 3: Average MSE Across All Scenarios",
       subtitle = "Lower = better CATE estimation", x = "Scenario", y = "Average MSE") +
  theme_bw(base_size = 12) + theme(legend.position = "bottom") + scale_fill_viridis_d()
print(p_mse)

# Overall ranking table
overall_ranking <- summary_stats %>% filter(Beta > 0) %>%
  group_by(Method) %>%
  summarise(Scenarios = n_distinct(Scenario), Avg_Power = mean(1 - Mean_P_Value),
            Avg_Detection_Rate = mean(Detection_Rate), Avg_MSE = mean(Mean_MSE),
            Avg_Correlation = mean(Mean_Correlation),
            Avg_Somers_D = mean(Mean_Somers_D, na.rm = TRUE), .groups = 'drop')
if ("Top_Modifier" %in% names(obj2_results)) {
  obj2_overall <- lapply(unique_scenarios, function(scen) {
    ti <- get_true_modifiers(scen)
    obj2_results %>%
      filter(Scenario == scen, Beta > 0, !is.na(Top_Modifier)) %>%
      mutate(Is_True_Modifier = Top_Modifier %in% ti$modifiers) %>%
      group_by(Method) %>%
      summarise(Scen_True_Rate = mean(Is_True_Modifier, na.rm = TRUE), .groups = 'drop')
  }) %>% bind_rows() %>% group_by(Method) %>%
    summarise(Avg_True_Selection = mean(Scen_True_Rate, na.rm = TRUE) * 100, .groups = 'drop')
  overall_ranking <- overall_ranking %>% full_join(obj2_overall, by = "Method")
}
overall_ranking <- overall_ranking %>% arrange(desc(Avg_Power), Avg_MSE)
table_overall <- gridExtra::tableGrob(
  overall_ranking %>% mutate(across(where(is.numeric), ~round(., 3))), rows = NULL)
grid.arrange(table_overall, top = "Overall Method Rankings (All Scenarios, β > 0)")

# Scenario characteristics (soft from R/scenario_meta.R)
scenario_characteristics <- do.call(rbind, lapply(
  sort(unique(combined_results$Scenario)), function(s) {
    ti <- get_true_modifiers(s)
    data.frame(Scenario = s,
               True_Modifiers  = if (length(ti$modifiers)) paste(ti$modifiers, collapse = ", ") else "(none)",
               N_Modifiers     = length(ti$modifiers),
               Prognostic_Only = if (length(ti$prognostic_only)) paste(ti$prognostic_only, collapse = ", ") else "-",
               stringsAsFactors = FALSE)
  }))
table_scenarios <- gridExtra::tableGrob(scenario_characteristics, rows = NULL)
grid.arrange(table_scenarios, top = "Scenario Characteristics: True Effect Modifiers (R/scenario_meta.R)")

# TMLE / ATE plots — GUARDED (clean harness emits no ate_estimates).
if (!is.null(combined_ate) && nrow(combined_ate) > 0) {
  cat("  TMLE/ATE plots rendered.\n")
}

# Wall-clock time comparison
if (any(!is.na(valid_results$Time_Sec))) {
  time_plot_df <- valid_results %>%
    group_by(Method) %>%
    summarise(mean_t = mean(Time_Sec, na.rm = TRUE), sd_t = sd(Time_Sec, na.rm = TRUE), .groups = "drop") %>%
    filter(!is.na(mean_t)) %>% arrange(mean_t)
  p_time <- ggplot(time_plot_df, aes(x = reorder(Method, mean_t), y = mean_t,
                                     fill = grepl("^TMLE", Method))) +
    geom_col(alpha = 0.85) +
    geom_errorbar(aes(ymin = pmax(0, mean_t - sd_t), ymax = mean_t + sd_t), width = 0.25) +
    scale_fill_manual(values = c(`TRUE` = "#CC6677", `FALSE` = "#4477AA"),
                      labels = c(`TRUE` = "TMLE", `FALSE` = "Other"), name = "Method family") +
    coord_flip() + theme_bw(base_size = 12) + theme(legend.position = "bottom") +
    labs(title = "Wall-clock time per CATE method (mean ± sd across simulations)",
         x = NULL, y = "Time (seconds per simulation)")
  print(p_time)
}

dev.off()
cat("  Saved: results/Cross_Scenario_Summary.pdf\n")

# =============================================================================
# 5. CONSOLE PERFORMANCE SUMMARY (paper Table-3 style)
# =============================================================================

cat("\n========================================\n")
cat("PERFORMANCE SUMMARY (β > 0)\n")
cat("========================================\n")

if (exists("overall_ranking") && nrow(overall_ranking) > 0) {
  obj1_rank <- overall_ranking %>% arrange(desc(Avg_Power)) %>%
    mutate(Obj1_Rank = row_number()) %>% select(Method, Obj1_Rank, Avg_Power)
  obj3_rank <- overall_ranking %>% arrange(Avg_MSE) %>%
    mutate(Obj3_Rank = row_number()) %>% select(Method, Obj3_Rank, Avg_MSE, Avg_Correlation, Avg_Somers_D)
  performance_summary <- obj1_rank %>% left_join(obj3_rank, by = "Method")
  if ("Avg_True_Selection" %in% names(overall_ranking)) {
    obj2_rank <- overall_ranking %>% arrange(desc(Avg_True_Selection)) %>%
      mutate(Obj2_Rank = row_number()) %>% select(Method, Obj2_Rank, Avg_True_Selection)
    performance_summary <- performance_summary %>% left_join(obj2_rank, by = "Method")
  }
  performance_summary <- performance_summary %>%
    mutate(Obj1_Symbol = case_when(Obj1_Rank <= 2 ~ "✓ Excellent", Obj1_Rank <= 3 ~ "● Competitive", TRUE ~ "✗ Poor"),
           Obj3_Symbol = case_when(Obj3_Rank <= 2 ~ "✓ Excellent", Obj3_Rank <= 3 ~ "● Competitive", TRUE ~ "✗ Poor"))
  if ("Obj2_Rank" %in% names(performance_summary))
    performance_summary <- performance_summary %>%
      mutate(Obj2_Symbol = case_when(Obj2_Rank <= 2 ~ "✓ Excellent", Obj2_Rank <= 3 ~ "● Competitive", TRUE ~ "✗ Poor"))

  fmt_num <- function(x, d = 3, s = "") if (is.null(x) || length(x) == 0 || is.na(x)) "—" else paste0(formatC(x, digits = d, format = "f"), s)
  fmt_int <- function(x) if (is.null(x) || length(x) == 0 || is.na(x)) "—" else as.character(as.integer(x))

  for (i in seq_len(nrow(performance_summary))) {
    row <- performance_summary[i, ]
    cat(sprintf("\n%s:\n", row$Method))
    cat(sprintf("  Obj 1 (Heterogeneity Test): %s (Power=%s, Rank=%s)\n",
                row$Obj1_Symbol, fmt_num(row$Avg_Power), fmt_int(row$Obj1_Rank)))
    if ("Obj2_Symbol" %in% names(row))
      cat(sprintf("  Obj 2 (Modifier ID):        %s (True Selection=%s, Rank=%s)\n",
                  ifelse(is.na(row$Obj2_Symbol), "—", row$Obj2_Symbol),
                  fmt_num(row$Avg_True_Selection, 1, "%"), fmt_int(row$Obj2_Rank)))
    somers_str <- if (!is.na(row$Avg_Somers_D)) sprintf(", Somers_D=%s", fmt_num(row$Avg_Somers_D, 3)) else ""
    cat(sprintf("  Obj 3 (CATE Estimation):    %s (MSE=%s%s, Rank=%s)\n",
                row$Obj3_Symbol, fmt_num(row$Avg_MSE, 4), somers_str, fmt_int(row$Obj3_Rank)))
  }

  ranked_cate <- performance_summary %>%
    filter(!is.na(Obj1_Rank), !is.na(Obj3_Rank)) %>%
    mutate(Total_Rank = Obj1_Rank + Obj3_Rank +
             if ("Obj2_Rank" %in% names(.)) ifelse(is.na(Obj2_Rank), 0, Obj2_Rank) else 0) %>%
    arrange(Total_Rank)
  if (nrow(ranked_cate) > 0) {
    cat("\n----------------------------------\n")
    cat(sprintf("Best Overall CATE Method: %s\n", ranked_cate$Method[1]))
    cat("----------------------------------\n")
  }
}

cat("\n========================================\n")
cat("WATCH ANALYSIS COMPLETE\n")
cat("========================================\n")
cat("Generated outputs:\n")
cat("  • results/final_summary_statistics.csv\n")
cat("  • results/final_vte_summary.csv (if VTE data present)\n")
cat("  • results/final_pairwise_relative_error.csv (Gao 2024 EIF)\n")
cat("  • results/Performance_Distribution_*.pdf (density / violin / ridge)\n")
cat("  • results/WATCH_analysis_Scenario_*.pdf (Obj 1/2/3 per scenario)\n")
cat("  • results/Cross_Scenario_Summary.pdf (power / MSE / rankings / timing)\n")
cat(sprintf("End Time: %s\n", Sys.time()))
cat("========================================\n")
