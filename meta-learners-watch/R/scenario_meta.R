# ============================================================================
# scenario_meta.R — SINGLE SOFT-CODED SOURCE OF TRUTH for scenario metadata
# (true effect modifiers, purely-prognostic vars, labels).
#
# Everything here is DERIVED from the actual data-generating definitions, so the
# analysis can never drift out of sync with the simulation:
#   - Scenarios 1-4  -> parsed from benchtm::scen_param's `pred`/`prog` strings
#                       (the exact Sun et al. 2022 DGP the pipeline samples from).
#   - Scenarios 5-8  -> parsed from the legacy custom definitions (mirror of
#                       R/data_generation.R::get_scenario_definition; only the
#                       legacy generate_X_dist path uses these).
#
# Sourced by scripts/aggregate.R and scripts/analyze_focus_methods.R. Only needs
# benchtm's scen_param dataset (no SuperLearner / model fitting).
#
# If you change which benchtm scenario maps to a number, or edit a legacy
# definition, the true modifiers / labels below update automatically — no
# hand-maintained lists to keep in sync.
# ============================================================================

# Pull every covariate name (X<number>) out of an expression string.
.extract_X_vars <- function(expr_string) {
  if (length(expr_string) == 0 || is.na(expr_string) || !nzchar(expr_string)) {
    return(character(0))
  }
  hits <- regmatches(expr_string, gregexpr("X[0-9]+", expr_string))[[1]]
  unique(hits)
}

# benchtm continuous scen_param (scenarios 1-4): one (prog, pred) per scenario.
.scen_meta_benchtm <- function() {
  e <- new.env()
  utils::data("scen_param", package = "benchtm", envir = e)
  sp <- get("scen_param", envir = e)
  sp <- sp[sp$type == "continuous", , drop = FALSE]
  sp$scenario <- rep(seq_len(nrow(sp) / 5L), each = 5L)   # 5 b1_rel levels / scenario
  sp[!duplicated(sp$scenario), c("scenario", "prog", "pred")]
}

# Legacy custom scenarios 5-8 — MIRRORS R/data_generation.R::get_scenario_definition.
# (Keep in sync if you edit the legacy DGP; the active benchtm scenarios 1-4 are
# parsed directly from scen_param and need no maintenance here.)
.scen_meta_legacy <- list(
  "5" = list(prog = "0.5*X3 + X7^2",                              pred = "X3^2 - 0.5*X3 + X7"),
  "6" = list(prog = "sin(2*pi*X3) + 0.5*X8",                      pred = "sin(pi*X3) * cos(pi*X7)"),
  "7" = list(prog = "exp(-2*abs(X3)) + X14",                      pred = "exp(-abs(X3-0.5))"),
  "8" = list(prog = "X3^2 + sin(2*pi*X7) - 0.5*log(abs(X14)+1)",  pred = "X3^3 - X3 + sin(pi*X7)")
)

# Return prog/pred strings for a scenario (benchtm for 1-4, legacy for 5-8).
.scen_expr <- function(scenario_num) {
  s <- as.numeric(scenario_num)
  if (s %in% 1:4) {
    sp  <- .scen_meta_benchtm()
    row <- sp[sp$scenario == s, , drop = FALSE]
    if (nrow(row) == 1) return(list(prog = row$prog, pred = row$pred))
  }
  key <- as.character(s)
  if (key %in% names(.scen_meta_legacy)) return(.scen_meta_legacy[[key]])
  list(prog = NA_character_, pred = NA_character_)
}

# PUBLIC: true effect modifiers (= vars in the predictive function) and
# purely-prognostic vars (in prog but not predictive), derived from the DGP.
get_true_modifiers <- function(scenario_num) {
  ex   <- .scen_expr(scenario_num)
  mods <- .extract_X_vars(ex$pred)
  list(modifiers       = mods,
       prognostic_only = setdiff(.extract_X_vars(ex$prog), mods))
}

# PUBLIC: soft scenario label, e.g. "S1 (X11)" or "S3 (X14, X1)".
scenario_label <- function(scenario_num) {
  mods <- get_true_modifiers(scenario_num)$modifiers
  sprintf("S%s (%s)", scenario_num,
          if (length(mods)) paste(mods, collapse = ", ") else "none")
}

# ============================================================================
# SHARED ANALYSIS CONFIG (single source for scripts/aggregate.R, analyze_focus_methods.R
# and plot_paper_figures.R — edit here once and every analysis honours it).
# ============================================================================

# Methods to drop from ALL analyses (e.g. EIF reports the raw pseudo-outcome).
EXCLUDE_METHODS <- c("DR-learner (Semiparametric EIF)")

# Consistent colour palette across every plot. Extra entries are harmless;
# methods not present in the data are simply unused.
METHOD_PAL <- c(
  "DINA Original"                   = "#1B9E77",
  "DINA Cross-Fit"                  = "#66C2A5",
  "DINA Gao (GLM k2)"               = "#00441B",
  "DINA Gao (LASSO k2)"             = "#41AB5D",
  "DINA Centered"                   = "#238B45",
  "DINA Centered Expanded"          = "#005A32",
  "DR-learner (Kennedy)"            = "#D95F02",
  "DR-learner (2-Model 2SL)"        = "#386CB0",
  "DR-learner (Semiparametric EIF)" = "#A6CEE3",
  "R-learner (LASSO)"               = "#E6AB02",
  "R-learner (Boost)"               = "#FFB300",
  "Dandl Extension"                 = "#7570B3",
  "WATCH (Residual LASSO)"          = "#E7298A",
  "WATCH (Residual Risk)"           = "#A6761D",
  "TMLE (Li 2026 VTE)"              = "#444444",
  "TMLE (Within-Trial Prog)"        = "#999999"
)
