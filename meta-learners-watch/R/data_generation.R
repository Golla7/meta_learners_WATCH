############################################################################
# R/data_generation.R
#
# Simulation data generation for the meta-learner comparison.
#
# Scenarios 1-4 reproduce the Sun et al. (2022) benchtm benchmark EXACTLY
# (Biometrical Journal, DOI 10.1002/bimj.202100337): synthpop covariates
# mimicking real Phase III inflammatory-disease trial data (p = 30; 8
# categorical, 22 numeric), an RCT treatment assignment, and a continuous
# outcome whose prognostic/predictive components and b0/b1 coefficients are
# pre-calibrated in benchtm's scen_param table. Scenarios 5-8 are legacy
# custom non-linear scenarios on parametric covariates, kept for optional
# extension runs; only those use the beta* calibration file.
#
# --- WORKFLOW (legacy scenarios 5-8 only) ---
# 1. CALIBRATE once per treatment scenario:
#    > rct_betas <- calculate_beta_star_for_scenarios(treatment_scenario = "RCT")
# 2. SAVE: > saveRDS(rct_betas, "data/beta_star_RCT.rds")
# 3. RUN:  scripts/run_parallel.sh validates the calibration file exists.
# For the benchtm scenarios 1-4 no calibration step is needed — b0/b1 come
# pre-calibrated from benchtm::scen_param.
############################################################################

# =====================================================
# 0. BENCHTM PACKAGE REQUIREMENT
# =====================================================
if (!require("benchtm", quietly = TRUE)) {
  stop(paste0(
    "\n=====================================\n",
    "ERROR: benchtm package is required!\n",
    "=====================================\n",
    "Install with: devtools::install_github('Sophie-Sun/benchtm')\n",
    "Then restart R and try again.\n"
  ))
}


# =====================================================
# 1. SCENARIO DEFINITIONS FOR THE LEGACY get_b() CALIBRATION PATH
# =====================================================

get_scenario_definition <- function(scenario_num) {
  scenario_num <- as.character(scenario_num)
  definitions <- list(
    "1" = list(prog = "0.5*(X3+X7)",      effect = "X3", b0 = 0),
    "2" = list(prog = "X14 - X8",         effect = "X14", b0 = 0.5),
    "3" = list(prog = "X3 - 0.5*X17",     effect = "as.numeric(X14 > 0.25 & X3 > 0.5)", b0 = 0),
    "4" = list(prog = "X11 - X14",        effect = "as.numeric(X14 > 0.3 | X4 > 0.5)", b0 = 0.2),
    "5" = list(prog = "0.5*X3 + X7^2",    effect = "X3^2 - 0.5*X3 + X7", b0 = 0),
    "6" = list(prog = "sin(2*pi*X3) + 0.5*X8", effect = "sin(pi*X3) * cos(pi*X7)", b0 = 0.2),
    "7" = list(prog = "exp(-2*abs(X3)) + X14", effect = "exp(-abs(X3-0.5))", b0 = 0.5),
    "8" = list(prog = "X3^2 + sin(2*pi*X7) - 0.5*log(abs(X14)+1)", effect = "X3^3 - X3 + sin(pi*X7)", b0 = 0)
  )
  if (!scenario_num %in% names(definitions)) stop("Invalid scenario number specified.")
  return(definitions[[scenario_num]])
}


# ============================================================================
# SUN ET AL. (2022) benchtm SCENARIOS (Biometrical Journal, DOI 10.1002/bimj.202100337)
# ----------------------------------------------------------------------------
# Scenarios 1-4 reproduce Table 1 / Table A.1 of the paper EXACTLY via benchtm:
#   - Covariates: generate_X_syn(n) -> p = 30 synthpop biomarkers mimicking real
#     Phase III inflammatory-disease data (8 categorical incl. X1/X4/X8 with
#     'Y'/'N' levels; 22 numeric scaled to [0,1]).
#   - The scenario forms AND the scaling factor s are baked into the prog/pred
#     strings in benchtm's scen_param table; b0/b1 are PRE-CALIBRATED there for
#     every beta1/beta* level (b1_rel in {0, 0.5, 1, 1.5, 2}) to give ~80%
#     interaction-test power and ~50% overall-effect power at n = 500. No
#     separate get_b() calibration step is required.
#   - beta_level in this pipeline == the paper's b1_rel. b0 changes with the
#     level (overall ATE held ~constant while heterogeneity grows), so we read
#     the whole (prog, pred, b0, b1) tuple straight from scen_param.
# ============================================================================

# Continuous rows of scen_param (scenarios 1-4, in order, 5 b1_rel levels each).
.benchtm_scen_continuous <- function() {
  envir <- new.env()
  utils::data("scen_param", package = "benchtm", envir = envir)
  sp <- get("scen_param", envir = envir)
  sp <- sp[sp$type == "continuous", , drop = FALSE]
  sp$scenario <- rep(1:4, each = 5)   # continuous block is laid out in scenario order
  sp
}

# Look up the (prog, pred, b0, b1) tuple for paper scenario (1-4) at a given
# beta_level (= b1_rel in {0, 0.5, 1, 1.5, 2}).
get_benchtm_scenario <- function(scenario, beta_level) {
  sp  <- .benchtm_scen_continuous()
  row <- sp[sp$scenario == scenario & abs(sp$b1_rel - beta_level) < 1e-8, , drop = FALSE]
  if (nrow(row) != 1L) {
    stop(sprintf(
      "No unique benchtm scenario for scenario=%s, beta_level(b1_rel)=%s. Available b1_rel: %s",
      scenario, beta_level, paste(sort(unique(sp$b1_rel)), collapse = ", ")
    ))
  }
  list(prog = row$prog, pred = row$pred, b0 = row$b0, b1 = row$b1, b1_rel = row$b1_rel)
}

# Generate one Sun et al. (2022) dataset: synthpop covariates + RCT treatment +
# continuous outcome with the pre-calibrated prog/pred/b0/b1. Returns the
# benchtm::generate_y data.frame (Y, trt, X1..X30, trt_effect = true CATE).
generate_paper_data <- function(n = 500, prog, pred, b0, b1, seed = 123,
                                treatment_scenario = "RCT") {
  set.seed(seed)
  if (!identical(treatment_scenario, "RCT")) {
    stop("Sun et al. (2022) scenarios (1-4) are defined for an RCT (P(A=1)=0.5) only.")
  }
  X   <- suppressWarnings(benchtm::generate_X_syn(n))
  trt <- benchtm::generate_trt(n, p_trt = 0.5, type = "exact")
  benchtm::generate_y(
    X, trt, prog = prog, pred = pred, b0 = b0, b1 = b1,
    type = "continuous", include_truth = TRUE, sigma_error = 1
  )
}


# ============================================================================
# CALIBRATION USING benchtm::get_b()
# ============================================================================

calculate_beta_star_for_scenarios <- function(scenarios_to_calc = 1:8,
                                              target_power = 0.8,
                                              alpha = 0.1,
                                              treatment_scenario = "RCT",
                                              n = 2000, 
                                              p = 20,
                                              scal = 100) {
  
  cat("========================================\n")
  cat("BETA* CALIBRATION USING benchtm::get_b()\n")
  cat("========================================\n\n")
  
  results_list <- list()
  
  for (i in seq_along(scenarios_to_calc)) {
    scenario <- scenarios_to_calc[i]
    cat(sprintf("--- Scenario %d (%s) ---\n", scenario, treatment_scenario))
    
    scenario_def <- get_scenario_definition(scenario)
    
    tryCatch({
      # Generate large covariate matrix (super-population)
      set.seed(scenario * 1000)
      X <- benchtm::generate_X_dist(n = n * scal, p = p, rho = 0.5)
      
      # Generate treatment based on scenario type
      if (treatment_scenario == "RCT") {
        trt <- benchtm::generate_trt(n = n * scal, p_trt = 0.5, type = "exact")
      } else if (treatment_scenario == "Observational_Simple") {
        trt <- benchtm::generate_trt(n = n * scal, type = "random", 
                                     X = X, prop = "0.5 * X5 - 0.5 * X10")
      } else if (treatment_scenario == "Observational_Complex") {
        trt <- benchtm::generate_trt(n = n * scal, type = "random", 
                                     X = X, 
                                     prop = "0.4*I(X5 > 0) - 0.4*I(X10 < 0) + 0.5 * X5 * X10")
      }
      
      # Call get_b() - calculates both b0 and b1
      b_params <- benchtm::get_b(
        X = X,
        scal = scal,
        prog = scenario_def$prog,
        pred = scenario_def$effect,
        trt = trt,
        type = "continuous",
        power = c(0.9, target_power),  # [overall test, interaction test]
        alpha = c(0.025, alpha),       # [overall test, interaction test]
        start = c(0, 0),
        sign_better = 1,
        sigma_error = 1,
        optim_method = "Nelder-Mead"
      )
      
      b1_calibrated <- b_params[2]  # This is your Beta_Star
      power_achieved <- attr(b_params, "power_results")
      
      cat(sprintf("  ✓ Beta_Star (b1) = %.4f\n", b1_calibrated))
      if (!is.null(power_achieved)) {
        cat(sprintf("  ✓ Power: Overall=%.3f, Interaction=%.3f\n",
                    power_achieved[1], power_achieved[2]))
      }
      
      results_list[[i]] <- data.frame(
        Scenario = scenario,
        Beta_Star = b1_calibrated
      )
      
    }, error = function(e) {
      cat(sprintf("  ✗ Error: %s\n", e$message))
      cat("  Using default value\n")
      
      defaults <- c(2.0, 1.5, 1.8, 2.5, 1.5, 2.0, 1.8, 2.0)
      results_list[[i]] <- data.frame(
        Scenario = scenario,
        Beta_Star = defaults[scenario]
      )
    })
    cat("\n")
  }
  
  results_df <- do.call(rbind, results_list)
  
  cat("--- FINAL CALIBRATION RESULTS ---\n")
  print(results_df)
  
  return(invisible(results_df))
}

# =====================================================
# 3. DATA GENERATION (Main function for simulation)
# =====================================================
DEFAULT_BETA_STAR_VALUES <- c(2.0, 1.5, 1.8, 2.5, 1.5, 2.0, 1.8, 2.0)

generate_simulation_data <- function(n = 500, p = 30, scenario = 1, beta_level = 1, seed = 123,
                                     beta_star_values = NULL,
                                     treatment_scenario = "RCT") {
  set.seed(seed)

  # --- Sun et al. (2022) scenarios (1-4): exact benchtm DGP. beta_level is the
  #     paper's b1_rel; (prog, pred, b0, b1) come PRE-CALIBRATED from scen_param,
  #     so beta_star_values / the calibration file are NOT used here. p is fixed
  #     at 30 by generate_X_syn (the p argument is ignored for these scenarios).
  if (scenario %in% 1:4) {
    sc <- get_benchtm_scenario(scenario, beta_level)
    return(generate_paper_data(n = n, prog = sc$prog, pred = sc$pred,
                               b0 = sc$b0, b1 = sc$b1, seed = seed,
                               treatment_scenario = treatment_scenario))
  }

  # --- Legacy custom scenarios (5-8): parametric generate_X_dist covariates with
  #     get_b()-calibrated beta* (b1 = beta_level * beta_star). ---
  scenario_def <- get_scenario_definition(scenario)

  beta_star_multiplier <- DEFAULT_BETA_STAR_VALUES[scenario]

  if (is.null(beta_star_values)) {
    warning("No calibration file provided. Using default beta* values.")
  } else if (scenario %in% beta_star_values$Scenario) {
    multiplier <- beta_star_values$Beta_Star[beta_star_values$Scenario == scenario]
    if (is.na(multiplier)) {
      warning(sprintf(
        "Scenario %d: Calibrated beta* is NA. Using default value %.2f",
        scenario, beta_star_multiplier
      ))
    } else if (multiplier <= 0) {
      warning(sprintf(
        "Scenario %d: Calibrated beta* is invalid (%.4f <= 0). Using default value %.2f",
        scenario, multiplier, beta_star_multiplier
      ))
    } else {
      if (multiplier > 10) {
        warning(sprintf(
          "Scenario %d: Calibrated beta* seems unusually large (%.4f). Using it anyway, but verify calibration.",
          scenario, multiplier
        ))
      }
      beta_star_multiplier <- multiplier
    }
  } else {
    warning(sprintf("Scenario %d not found in provided calibration file. Using default.", scenario))
  }

  b1 <- beta_level * beta_star_multiplier

  data <- generate_comparison_data(n = n, p = p, seed = seed,
                                   prog_function = scenario_def$prog,
                                   effect_modifier = scenario_def$effect,
                                   b0 = scenario_def$b0, b1 = b1,
                                   treatment_scenario = treatment_scenario)
  return(data)
}

# =====================================================
# 3.1 CORE DATA GENERATION WRAPPER (USES BENCHTM)
# =====================================================
generate_comparison_data <- function(n = 2000, p = 20, seed = 123,
                                     prog_function, effect_modifier,
                                     b0, b1,
                                     sigma_error = 1, rho = 0.5,
                                     treatment_scenario = c("RCT", "Observational_Simple", "Observational_Complex")) {

  treatment_scenario <- match.arg(treatment_scenario)
  set.seed(seed)

  X_numeric <- benchtm::generate_X_dist(n, p, rho = rho)

  if (treatment_scenario == "RCT") {
    trt <- benchtm::generate_trt(n, p_trt = 0.5, type = "exact")
  } else if (treatment_scenario == "Observational_Simple") {
    prop_formula <- "0.5 * X5 - 0.5 * X10"
    trt <- benchtm::generate_trt(n, type = "random", X = X_numeric, prop = prop_formula)
  } else if (treatment_scenario == "Observational_Complex") {
    prop_formula <- "0.4*I(X5 > 0) - 0.4*I(X10 < 0) + 0.5 * X5 * X10"
    trt <- benchtm::generate_trt(n, type = "random", X = X_numeric, prop = prop_formula)
  }

  final_data <- benchtm::generate_y(
    X_numeric,
    trt,
    prog = prog_function,
    pred = effect_modifier,
    b0 = b0,
    b1 = b1,
    type = "continuous",
    include_truth = TRUE,
    sigma_error = sigma_error
  )

  return(final_data)
}
