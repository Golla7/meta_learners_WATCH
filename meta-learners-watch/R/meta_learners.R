############################################################################
# R/meta_learners.R
#
# A SAME-FOOTING comparison harness for HTE meta-learners, built around the
# three objectives of the WATCH workflow (Sechidis et al., arXiv:2502.00713):
#   Objective 1 — global test for treatment-effect heterogeneity
#   Objective 2 — ranking covariates by effect-modification strength
#   Objective 3 — estimating the individualized treatment effect (CATE)
#
# DESIGN PRINCIPLE: "SAME FOOTING"
# --------------------------------
# When meta-learners each fit their own stage-1 nuisance models, "which
# estimator is better" is confounded with "which got better nuisances".
# Here ONE unified stage-1 estimator (estimate_nuisances, §2) reproduces the
# reference DR-learner's nuisance configuration (SL.glmnet + SL.cforest /
# method.CC_LS outcome models; SL.glmnet / method.NNLS2 propensity; 5-fold
# cross-fitting stratified on Y; propensity estimated even in an RCT), and
# EVERY meta-learner consumes those SAME nuisances for its second stage.
# Methods therefore differ ONLY in the stage-2 CATE regression, so
# comparisons are clean.
#
# Reference learner ...... DR-learner (2-Model 2SL), after K. Sechidis
# Comparison methods ..... DINA Original, DINA Cross-Fit, DINA Centered,
#                          DINA Centered Expanded, DR-learner (Kennedy),
#                          R-learner (LASSO), Dandl Extension,
#                          WATCH (Residual LASSO), WATCH (Residual Risk),
#                          TMLE (Li 2026 VTE)
# Self-contained pair .... DINA Gao (GLM k2) + DINA Gao (LASSO k2) — faithful
#                          Gao & Hastie (2022) replicas. Algorithm 1 couples
#                          the nuisance split to the causal split, so each
#                          runs its own k=2 cross-fit (same SuperLearner
#                          library as estimate_nuisances; marginal
#                          nu(x)=E[Y|X] per Gao's Gaussian shortcut); they
#                          differ only in stage 2 — UNPENALIZED lm (paper
#                          §4.2; == linear R-learner for Gaussian Y) vs L1
#                          cv.glmnet (paper footnote). These two are the ONLY
#                          learners not on the shared nuisances.
#
# Every learner returns the same standardized list (tau, tau_eval, global-test
# p-values, variable importance), and compare_meta_learners() (§4) dispatches
# all of them on one shared nuisance object. The guarded demo in §5 runs when
# this file is executed directly from the repo root:
#   Rscript R/meta_learners.R
############################################################################


# ==========================================================================
# §0. LIBRARIES & PACKAGE BOOTSTRAP
# ==========================================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(SuperLearner)
  library(nnls)
  library(caret)
  library(glmnet)
  library(party)        # SL.cforest base learner + cforest VI
  library(partykit)     # ctree_control for the Dandl pmforest
  library(model4you)    # pmforest / pmodel (Dandl)
  library(coin)         # global heterogeneity test
  library(permimp)      # permutation variable importance
  library(sandwich)     # estfun() for the WATCH score residual
  library(Hmisc)        # rcorr.cens() for Somers' D
})

# rlearner (R-learner LASSO) is installed from GitHub.
if (!requireNamespace("rlearner", quietly = TRUE)) {
  message("meta_learners: 'rlearner' not installed; ",
          "R-learner (LASSO) will be skipped. ",
          "Install via devtools::install_github('xnie/rlearner').")
}

# TE-VIM TMLE (Li 2026, Psi_2 / Psi_3) — faithful port of HaodongL/te_vim.
# Defines compute_te_vim() + TMLE_VIM/TMLE_VTE/compute_psi3. Sourced relative
# to the repo root (all scripts in this project run from the repo root).
if (file.exists("R/te_vim.R")) {
  source("R/te_vim.R")
} else if (file.exists("te_vim.R")) {
  source("te_vim.R")   # fallback when sourcing from inside R/
} else {
  message("meta_learners: 'R/te_vim.R' not found; ",
          "compute_vim=TRUE will be unavailable.")
}


# ==========================================================================
# §1. SHARED UTILITIES
#
# The small helpers every learner needs, so this file runs standalone.
# ==========================================================================

# --------------------------------------------------------------------------
# global_heterogeneity_test() — Sechidis et al. (arXiv:2502.00713), Objective 1.
#
# Tests H0: P(psi | X) = P(psi) via coin::independence_test, returning both the
# max-type and quadratic-type permutation p-values. `pseudo_outcome` is the
# per-observation individual-treatment-effect proxy appropriate to the caller.
# --------------------------------------------------------------------------
global_heterogeneity_test <- function(pseudo_outcome, X) {
  na_result <- list(p_max = NA_real_, p_quad = NA_real_)
  if (is.null(pseudo_outcome) || length(pseudo_outcome) < 10) return(na_result)

  pseudo_outcome <- as.numeric(pseudo_outcome)
  Xd <- as.data.frame(X)
  if (nrow(Xd) != length(pseudo_outcome)) return(na_result)

  valid <- is.finite(pseudo_outcome)
  if (sum(valid) < 10) return(na_result)

  y_v <- pseudo_outcome[valid]
  X_v <- Xd[valid, , drop = FALSE]

  for (col in names(X_v)) {
    if (is.character(X_v[[col]])) X_v[[col]] <- as.factor(X_v[[col]])
  }
  # Drop zero-variance columns (coin errors on them; they carry no signal).
  keep <- vapply(X_v, function(col) {
    if (is.numeric(col)) sd(col, na.rm = TRUE) > 1e-12
    else length(unique(col[!is.na(col)])) > 1L
  }, logical(1))
  if (!any(keep)) return(na_result)
  X_v <- X_v[, keep, drop = FALSE]
  if (sd(y_v, na.rm = TRUE) < 1e-12) return(na_result)

  test_df <- data.frame(.gh_y = y_v, X_v, check.names = FALSE)
  rhs     <- paste(sprintf("`%s`", setdiff(names(test_df), ".gh_y")), collapse = " + ")
  formula <- as.formula(paste(".gh_y ~", rhs))

  p_max <- tryCatch({
    t_obj <- coin::independence_test(formula, data = test_df, teststat = "maximum")
    suppressWarnings(as.numeric(coin::pvalue(t_obj)))
  }, error = function(e) NA_real_)
  p_quad <- tryCatch({
    t_obj <- coin::independence_test(formula, data = test_df, teststat = "quadratic")
    suppressWarnings(as.numeric(coin::pvalue(t_obj)))
  }, error = function(e) NA_real_)

  list(p_max = p_max, p_quad = p_quad)
}

# --------------------------------------------------------------------------
# compute_somers_d()
# Somers' Dxy rank correlation between estimated and true CATE — the
# Objective-3 ranking metric recommended by Sechidis et al. We use
# Hmisc::rcorr.cens (not somers2, which only handles binary outcomes).
# --------------------------------------------------------------------------
compute_somers_d <- function(estimated, true) {
  if (length(estimated) < 2 || length(estimated) != length(true)) return(NA_real_)
  valid <- !is.na(estimated) & !is.na(true)
  if (sum(valid) < 2) return(NA_real_)
  est_v <- estimated[valid]; true_v <- true[valid]
  if (sd(est_v) < 1e-12 || sd(true_v) < 1e-12) return(NA_real_)

  out <- tryCatch(
    Hmisc::rcorr.cens(est_v, true_v),
    error   = function(e) NULL,
    warning = function(w) suppressWarnings(Hmisc::rcorr.cens(est_v, true_v))
  )
  if (!is.null(out) && "Dxy" %in% names(out) && is.finite(out[["Dxy"]])) {
    return(unname(out[["Dxy"]]))
  }
  tryCatch(unname(cor(est_v, true_v, method = "kendall")), error = function(e) NA_real_)
}

# --------------------------------------------------------------------------
# evaluate_method()
# CATE accuracy metrics vs. the oracle tau_true.
# --------------------------------------------------------------------------
evaluate_method <- function(tau_hat, tau_true) {
  valid_idx <- !is.na(tau_hat) & !is.na(tau_true)
  tau_hat <- tau_hat[valid_idx]; tau_true <- tau_true[valid_idx]
  if (length(tau_hat) == 0) {
    return(list(MSE = NA, Bias = NA, Variance = NA, Correlation = NA,
                MAE = NA, RMSE = NA, Somers_D = NA))
  }
  mse <- mean((tau_hat - tau_true)^2); bias <- mean(tau_hat - tau_true)
  variance <- var(tau_hat - tau_true)
  correlation <- if (sd(tau_hat) < 1e-12) NA_real_ else cor(tau_hat, tau_true)
  list(MSE = mse, Bias = bias, Variance = variance, Correlation = correlation,
       MAE = mean(abs(tau_hat - tau_true)), RMSE = sqrt(mse),
       Somers_D = compute_somers_d(tau_hat, tau_true))
}

# --------------------------------------------------------------------------
# .dina_coef_vi()
# Coefficient variable importance = |beta_j| * sd(x_j), random tie-break.
# --------------------------------------------------------------------------
.dina_coef_vi <- function(beta_final, X_matrix) {
  modifier_cols <- setdiff(colnames(X_matrix), "(Intercept)")
  if (length(modifier_cols) == 0) {
    return(data.frame(variable = character(0), importance = numeric(0),
                      stringsAsFactors = FALSE))
  }
  col_sd <- apply(X_matrix[, modifier_cols, drop = FALSE], 2, stats::sd)
  imp <- abs(beta_final[modifier_cols]) * col_sd
  imp[!is.finite(imp)] <- 0
  vi <- data.frame(variable = modifier_cols, importance = as.numeric(imp),
                   stringsAsFactors = FALSE)
  vi[order(-vi$importance, sample(nrow(vi))), ]  # random tie-break
}

# --------------------------------------------------------------------------
# compute_relative_error_eif() — Gao (2024) pairwise relative-error EIF.
#
# VERIFIED EXACT MATCH to the official repo
# ZijunGao/Causal-validation-IF (manuscript/code for submission/helper.R):
#   efficient.relative.MSE() = mean( tau1^2 - tau2^2 - 2*weighted.ATE(...) )
#   weighted.ATE(weight, identity link) =
#       weight * ( W*(Y-mu1)/pi + mu1 - (1-W)*(Y-mu0)/(1-pi) - mu0 )
# Expanding with weight = (tau1 - tau2) reproduces the `psi` below exactly.
# The only addition here is a defensive clip of the propensity to [0.01,0.99].
#
# Y, W, e_hat/mu0_hat/mu1_hat MUST come from the held-out EVALUATION sample
# (out-of-sample nuisances), and tau1/tau2 are the two estimators' predictions
# on that same eval sample.
# --------------------------------------------------------------------------
compute_relative_error_eif <- function(tau1, tau2, Y, W,
                                        e_hat, mu0_hat, mu1_hat, alpha = 0.05) {
  na_result <- list(delta = NA_real_, se = NA_real_, ci_lo = NA_real_,
                    ci_hi = NA_real_, ci_excludes_zero = NA, n = 0L)
  if (length(tau1) != length(tau2) || length(tau1) != length(Y)) return(na_result)
  valid <- !is.na(tau1) & !is.na(tau2) & !is.na(Y) & !is.na(W) &
           !is.na(e_hat) & !is.na(mu0_hat) & !is.na(mu1_hat) &
           is.finite(tau1) & is.finite(tau2) & is.finite(e_hat)
  if (sum(valid) < 10) { na_result$n <- sum(valid); return(na_result) }

  t1 <- as.numeric(tau1[valid]); t2 <- as.numeric(tau2[valid])
  Yv <- as.numeric(Y[valid]);    Wv <- as.numeric(W[valid])
  e  <- pmax(0.01, pmin(as.numeric(e_hat[valid]), 0.99))
  m0 <- as.numeric(mu0_hat[valid]); m1 <- as.numeric(mu1_hat[valid])
  n  <- length(t1)

  # Equation (7) per-observation EIF (the "- delta" centering cancels in sd).
  pseudo_dr <- (Wv * (Yv - m1)) / e + m1 - ((1 - Wv) * (Yv - m0)) / (1 - e) - m0
  psi <- (t1^2 - t2^2) - 2 * (t1 - t2) * pseudo_dr

  delta <- mean(psi); se <- sd(psi) / sqrt(n)
  if (!is.finite(se) || se == 0) {
    return(list(delta = delta, se = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
                ci_excludes_zero = NA, n = n))
  }
  z <- qnorm(1 - alpha / 2)
  list(delta = delta, se = se, ci_lo = delta - z * se, ci_hi = delta + z * se,
       ci_excludes_zero = (delta - z * se > 0) || (delta + z * se < 0), n = n)
}

# Small internal helpers.
.coerce_numeric_df <- function(X) {
  X <- as.data.frame(X)
  for (col in names(X)) if (!is.numeric(X[[col]])) X[[col]] <- as.numeric(as.factor(X[[col]]))
  X
}
.r_learner_pseudo <- function(Y, trt, e_hat, mu0_hat, mu1_hat, thresh = 0.01) {
  if (is.null(e_hat) || is.null(mu0_hat) || is.null(mu1_hat)) return(NULL)
  m_hat   <- e_hat * mu1_hat + (1 - e_hat) * mu0_hat
  w_resid <- as.numeric(trt) - e_hat
  out     <- rep(NA_real_, length(Y))
  ok      <- abs(w_resid) > thresh
  out[ok] <- (as.numeric(Y)[ok] - m_hat[ok]) / w_resid[ok]
  out
}
.num_or_na <- function(x) if (is.null(x)) NA_real_ else suppressWarnings(as.numeric(x))


# ==========================================================================
# §2. UNIFIED NUISANCE ESTIMATOR  (the heart of "same footing")
#
# Reproduces the EXACT stage-1 config of the reference DR-learner
# (2-Model 2SL, after K. Sechidis):
#   - 5-fold cross-fitting, folds STRATIFIED ON Y (caret::createFolds).
#   - mu0, mu1 : per-arm SuperLearner(SL.glmnet, SL.cforest), method.CC_LS.
#   - pi       : SuperLearner(SL.glmnet), method.NNLS2 — ESTIMATED EVEN IN AN
#                RCT (the Kostas convention), then clipped to [0.01, 0.99].
#
# Every meta-learner in §3 consumes the SAME returned nuisances, so methods
# differ only in their stage-2 CATE regression.
#
# `light_config = TRUE` swaps SL.cforest -> SL.ranger for fast smoke tests.
# The shared DR (AIPW) pseudo-outcome `phi` is computed once and returned so
# DR-style learners (reference, Kennedy) need not recompute it.
# ==========================================================================
estimate_nuisances <- function(X, Y, trt, k = 5, seed = NULL,
                               treatment_scenario = "RCT",
                               light_config = FALSE, verbose = FALSE) {
  X   <- .coerce_numeric_df(X)
  Y   <- as.numeric(Y); trt <- as.numeric(trt)
  n   <- nrow(X)

  outcome_lib <- if (light_config) list("SL.glmnet", "SL.ranger")
                 else                list("SL.glmnet", "SL.cforest")
  prop_lib    <- "SL.glmnet"

  if (!is.null(seed)) set.seed(seed)
  # createFolds(list = FALSE) returns an integer fold id per row, binned on Y
  # (stratification on the outcome, matching the reference learner).
  fold_ids <- caret::createFolds(y = Y, k = k, list = FALSE)
  uf <- sort(unique(fold_ids))

  e_hat <- numeric(n); mu0_hat <- numeric(n); mu1_hat <- numeric(n)

  if (verbose) cat(sprintf("estimate_nuisances: %d-fold, %s outcomes, est. propensity\n",
                           length(uf), if (light_config) "glmnet+ranger" else "glmnet+cforest"))

  for (k_i in uf) {
    test_idx  <- which(fold_ids == k_i)
    train_idx <- which(fold_ids != k_i)
    tr_ctrl   <- intersect(train_idx, which(trt == 0))
    tr_trt    <- intersect(train_idx, which(trt == 1))

    # --- mu1: outcome model, treated arm (CC_LS) ---
    mu1_hat[test_idx] <- if (length(tr_trt) < 5) mean(Y[tr_trt]) else tryCatch({
      f <- SuperLearner::SuperLearner(
        Y = Y[tr_trt], X = X[tr_trt, , drop = FALSE],
        SL.library = outcome_lib, family = "gaussian",
        method = "method.CC_LS", verbose = FALSE)
      as.numeric(predict(f, newdata = X[test_idx, , drop = FALSE])$pred)
    }, error = function(e) rep(mean(Y[tr_trt]), length(test_idx)))

    # --- mu0: outcome model, control arm (CC_LS) ---
    mu0_hat[test_idx] <- if (length(tr_ctrl) < 5) mean(Y[tr_ctrl]) else tryCatch({
      f <- SuperLearner::SuperLearner(
        Y = Y[tr_ctrl], X = X[tr_ctrl, , drop = FALSE],
        SL.library = outcome_lib, family = "gaussian",
        method = "method.CC_LS", verbose = FALSE)
      as.numeric(predict(f, newdata = X[test_idx, , drop = FALSE])$pred)
    }, error = function(e) rep(mean(Y[tr_ctrl]), length(test_idx)))

    # --- pi: propensity (NNLS2). Estimated even for RCT (Kostas convention). ---
    e_hat[test_idx] <- tryCatch({
      f <- SuperLearner::SuperLearner(
        Y = trt[train_idx], X = X[train_idx, , drop = FALSE],
        SL.library = prop_lib, family = "binomial",
        method = "method.NNLS2", verbose = FALSE)
      as.numeric(predict(f, newdata = X[test_idx, , drop = FALSE])$pred)
    }, error = function(e) rep(mean(trt[train_idx]), length(test_idx)))
  }

  e_hat <- pmax(0.01, pmin(e_hat, 0.99))

  # Shared DR / AIPW pseudo-outcome, computed once from the cross-fitted
  # nuisances (identical to .dr_pseudo / the reference's per-fold phi).
  m_hat <- mu0_hat * (1 - trt) + mu1_hat * trt
  phi   <- mu1_hat - mu0_hat + (trt - e_hat) / (e_hat * (1 - e_hat)) * (Y - m_hat)

  list(e_hat = e_hat, mu0_hat = mu0_hat, mu1_hat = mu1_hat,
       fold_ids = fold_ids, phi = phi,
       config = list(k = length(uf), outcome_lib = outcome_lib,
                     prop_lib = prop_lib, light = light_config))
}

# --------------------------------------------------------------------------
# estimate_nuisances_eval() — companion for the held-out Gao-2024 EIF grid.
#
# X_eval is INDEPENDENT of training, so no fold splitting: fit once on ALL
# training data and predict at X_eval. Same SuperLearner config as
# estimate_nuisances(). Returns honest out-of-sample (e_hat, mu0_hat, mu1_hat)
# at the eval points for compute_relative_error_eif().
# --------------------------------------------------------------------------
estimate_nuisances_eval <- function(X, Y, trt, X_eval, seed = NULL,
                                    treatment_scenario = "RCT",
                                    light_config = FALSE) {
  X      <- .coerce_numeric_df(X)
  X_eval <- .coerce_numeric_df(X_eval)
  Y <- as.numeric(Y); trt <- as.numeric(trt); n_eval <- nrow(X_eval)

  outcome_lib <- if (light_config) list("SL.glmnet", "SL.ranger")
                 else                list("SL.glmnet", "SL.cforest")
  if (!is.null(seed)) set.seed(seed)

  ctrl_idx <- which(trt == 0); trt_idx <- which(trt == 1)

  e_hat <- tryCatch({
    f <- SuperLearner::SuperLearner(Y = trt, X = X, SL.library = "SL.glmnet",
                                    family = "binomial", method = "method.NNLS2",
                                    verbose = FALSE)
    as.numeric(predict(f, newdata = X_eval)$pred)
  }, error = function(e) rep(mean(trt), n_eval))
  e_hat <- pmax(0.01, pmin(e_hat, 0.99))

  mu1_hat <- if (length(trt_idx) < 5) rep(mean(Y[trt_idx]), n_eval) else tryCatch({
    f <- SuperLearner::SuperLearner(Y = Y[trt_idx], X = X[trt_idx, , drop = FALSE],
                                    SL.library = outcome_lib, family = "gaussian",
                                    method = "method.CC_LS", verbose = FALSE)
    as.numeric(predict(f, newdata = X_eval)$pred)
  }, error = function(e) rep(mean(Y[trt_idx]), n_eval))

  mu0_hat <- if (length(ctrl_idx) < 5) rep(mean(Y[ctrl_idx]), n_eval) else tryCatch({
    f <- SuperLearner::SuperLearner(Y = Y[ctrl_idx], X = X[ctrl_idx, , drop = FALSE],
                                    SL.library = outcome_lib, family = "gaussian",
                                    method = "method.CC_LS", verbose = FALSE)
    as.numeric(predict(f, newdata = X_eval)$pred)
  }, error = function(e) rep(mean(Y[ctrl_idx]), n_eval))

  list(e_hat = e_hat, mu0_hat = mu0_hat, mu1_hat = mu1_hat)
}


# ==========================================================================
# §3. META-LEARNERS (reference + focus-8)
#
# UNIFORM SIGNATURE:  learner_xxx(X, Y, trt, nuis, X_eval = NULL,
#                                 seed = NULL, nperm = 10)
#   - `nuis` is the §2 estimate_nuisances() output (shared by ALL learners).
#   - returns the STANDARDIZED list:
#       list(tau, tau_eval, method, p_value_max, p_value_quad, vi,
#            [vi_coef], [vte fields])
#
# DATA-LEAKAGE CONSTRAINT: learners receive ONLY X_eval (never Y_eval/trt_eval)
# — they use it solely to emit out-of-sample tau_eval predictions. Eval-set
# outcomes/treatment are used exclusively by the §4 driver's EIF step.
# ==========================================================================

# --------------------------------------------------------------------------
# 1. REFERENCE — DR-learner (2-Model 2SL), after K. Sechidis.
#
# Stage 1 now comes from the SHARED `nuis` (the private 5-fold refit is
# dropped — this is the one intended deviation from the byte-exact replica,
# made so every method shares stage-1). Stage 2 is kept faithful: for each of
# the K folds, fit SuperLearner(glmnet, cforest; CC_LS) on that fold's phi,
# predict on ALL N, and AVERAGE over folds. Global test on the out-of-fold
# phi; VI via cforest(ntree = 2000) + permimp.
# --------------------------------------------------------------------------
learner_dr_2model_2sl <- function(X, Y, trt, nuis, X_eval = NULL,
                                  seed = NULL, nperm = 20) {
  X  <- .coerce_numeric_df(X); N <- nrow(X)
  Xe <- if (!is.null(X_eval)) .coerce_numeric_df(X_eval) else NULL
  phi <- nuis$phi; fold_ids <- nuis$fold_ids
  uf  <- sort(unique(fold_ids))
  out_lib <- nuis$config$outcome_lib

  tau      <- numeric(N)
  tau_eval <- if (!is.null(Xe)) numeric(nrow(Xe)) else NULL
  n_used   <- 0L
  for (k_i in uf) {
    test_ids <- which(fold_ids == k_i)
    if (length(test_ids) < 3) next
    SL_phi <- tryCatch(SuperLearner::SuperLearner(
      Y = phi[test_ids], X = X[test_ids, , drop = FALSE],
      SL.library = out_lib, family = "gaussian",
      method = "method.CC_LS", verbose = FALSE), error = function(e) NULL)
    if (is.null(SL_phi)) next
    tau <- tau + as.numeric(predict(SL_phi, newdata = X)$pred)
    if (!is.null(Xe)) tau_eval <- tau_eval + as.numeric(predict(SL_phi, newdata = Xe)$pred)
    n_used <- n_used + 1L
  }
  if (n_used == 0L) { tau <- rep(mean(phi, na.rm = TRUE), N)
                      if (!is.null(Xe)) tau_eval <- rep(mean(phi, na.rm = TRUE), nrow(Xe)) }
  else { tau <- tau / n_used; if (!is.null(Xe)) tau_eval <- tau_eval / n_used }

  het <- global_heterogeneity_test(phi, X)
  vi  <- tryCatch({
    fit <- party::cforest(y ~ ., data.frame(X, y = phi),
                          control = party::cforest_unbiased(ntree = 2000))
    cf  <- permimp::permimp(fit, conditional = FALSE, nperm = nperm)
    df  <- data.frame(variable = names(cf$values), importance = cf$values,
                      stringsAsFactors = FALSE)
    df[order(-df$importance), ]
  }, error = function(e)
    data.frame(variable = names(X), importance = rep(0, ncol(X)), stringsAsFactors = FALSE))

  list(tau = tau, tau_eval = tau_eval, method = "DR-learner (2-Model 2SL)",
       p_value = het$p_max, p_value_max = het$p_max, p_value_quad = het$p_quad, vi = vi)
}

# --------------------------------------------------------------------------
# 2. DR-learner (Kennedy), after Kennedy (2023).
# Stage 2: K-fold cross-fit a SuperLearner on the shared phi (OOF predictions);
# a final full-data fit gives tau_eval. Stage-2 library unified to the
# reference's (glmnet, cforest) for footing parity.
# --------------------------------------------------------------------------
learner_dr_kennedy <- function(X, Y, trt, nuis, X_eval = NULL,
                               seed = NULL, nperm = 10) {
  X <- .coerce_numeric_df(X); n <- nrow(X)
  phi <- nuis$phi; fold_ids <- nuis$fold_ids
  uf  <- sort(unique(fold_ids)); out_lib <- nuis$config$outcome_lib
  tau_hat <- rep(NA_real_, n)

  for (k_i in uf) {
    train_idx <- which(fold_ids != k_i); test_idx <- which(fold_ids == k_i)
    if (length(train_idx) < 10 || length(test_idx) < 1) next
    fit <- tryCatch(SuperLearner::SuperLearner(
      Y = phi[train_idx], X = X[train_idx, , drop = FALSE],
      SL.library = out_lib, family = "gaussian", method = "method.CC_LS",
      verbose = FALSE), error = function(e) NULL)
    tau_hat[test_idx] <- if (is.null(fit)) mean(phi[train_idx], na.rm = TRUE) else
      tryCatch(as.numeric(predict(fit, X[test_idx, , drop = FALSE])$pred),
               error = function(e) rep(mean(phi[train_idx], na.rm = TRUE), length(test_idx)))
  }
  if (any(is.na(tau_hat))) tau_hat[is.na(tau_hat)] <- mean(phi, na.rm = TRUE)

  tau_eval <- NULL
  if (!is.null(X_eval)) {
    Xe <- .coerce_numeric_df(X_eval)
    full <- tryCatch(SuperLearner::SuperLearner(
      Y = phi, X = X, SL.library = out_lib, family = "gaussian",
      method = "method.CC_LS", verbose = FALSE), error = function(e) NULL)
    tau_eval <- if (is.null(full)) rep(mean(phi, na.rm = TRUE), nrow(Xe)) else
      tryCatch(as.numeric(predict(full, Xe)$pred),
               error = function(e) rep(mean(phi, na.rm = TRUE), nrow(Xe)))
  }

  het <- global_heterogeneity_test(phi, X)
  vi  <- tryCatch({
    cf_fit <- party::cforest(phi ~ ., data = data.frame(X, phi = phi),
                             control = party::cforest_unbiased(
                               mtry = ceiling(sqrt(ncol(X))), ntree = 200))
    v  <- permimp::permimp(cf_fit, conditional = FALSE, nperm = nperm)
    df <- data.frame(variable = names(v$values), importance = v$values, stringsAsFactors = FALSE)
    df[order(-df$importance), ]
  }, error = function(e)
    data.frame(variable = names(X), importance = rep(0, ncol(X)), stringsAsFactors = FALSE))

  list(tau = as.vector(tau_hat), tau_eval = tau_eval, method = "DR-learner (Kennedy)",
       p_value = het$p_max, p_value_max = het$p_max, p_value_quad = het$p_quad, vi = vi)
}

# --------------------------------------------------------------------------
# Shared DINA design + heterogeneity/VI block (used by DINA Original & Cross-Fit).
# --------------------------------------------------------------------------
# HARMONIZED ENCODING (2026-06-11): factors are coerced to single numeric
# columns via .coerce_numeric_df — the SAME encoding learner_rlearner_lasso
# uses — instead of model.matrix dummy expansion. X3 is the only 5-level
# factor in the benchtm DGP; dummy expansion gave it ~4 lasso columns (~4x
# the null selection probability of every other covariate; X3 top-modifier
# rate ~30% under beta=0), confounding Obj 2/3 with the encoding choice.
.dina_design <- function(X, W_centered) {
  X_matrix <- cbind(`(Intercept)` = 1, as.matrix(.coerce_numeric_df(X)))
  WX_data  <- as.data.frame(W_centered * X_matrix)
  colnames(WX_data) <- c("resW", paste0("resW.", colnames(X_matrix)[-1]))
  list(X_matrix = X_matrix, WX_matrix = as.matrix(WX_data),
       pen_factor = c(0, rep(1, ncol(X_matrix) - 1)))
}
.dina_eval_tau <- function(X_eval, X_matrix, beta_final, fallback) {
  Xe <- as.data.frame(X_eval)
  tryCatch({
    Xe_m  <- cbind(`(Intercept)` = 1, as.matrix(.coerce_numeric_df(Xe)))
    keep  <- intersect(colnames(X_matrix), colnames(Xe_m))
    pos   <- match(keep, colnames(X_matrix))
    as.vector(Xe_m[, keep, drop = FALSE] %*% beta_final[pos])
  }, error = function(e) rep(fallback, nrow(Xe)))
}
.dina_het_vi <- function(X, Y_residual, W_centered, tau_hat, nperm) {
  pseudo_tau <- rep(NA, length(W_centered)); nz <- abs(W_centered) > 0.01
  pseudo_tau[nz] <- Y_residual[nz] / W_centered[nz]
  het <- list(p_max = NA_real_, p_quad = NA_real_)
  ok  <- !is.na(pseudo_tau) & is.finite(pseudo_tau)
  if (sum(ok) > 10) het <- global_heterogeneity_test(pseudo_tau[ok], X[ok, , drop = FALSE])
  vi <- tryCatch({
    cf <- party::cforest(tau_hat ~ ., data = data.frame(X, tau_hat = tau_hat),
                         control = party::cforest_unbiased(
                           mtry = ceiling(sqrt(ncol(X))), ntree = 200))
    v  <- permimp::permimp(cf, conditional = FALSE, nperm = nperm)
    df <- data.frame(variable = names(v$values), importance = v$values, stringsAsFactors = FALSE)
    # Random tie-break (as in .dina_coef_vi). Without it, a ~constant tau_hat
    # gives all-zero importances and order() keeps column order — the first
    # covariate (X1) was crowned Top_Modifier deterministically under the null.
    df[order(-df$importance, sample(nrow(df))), ]
  }, error = function(e)
    data.frame(variable = names(X), importance = rep(0, ncol(X)), stringsAsFactors = FALSE))
  list(het = het, vi = vi)
}

# --------------------------------------------------------------------------
# 3. DINA Original — Gao & Hastie (2022), single full-sample stage-2 LASSO.
#
# !!! SAMPLE-SPLITTING CAVEAT (intentional) !!!
# This "naive" variant fits a SINGLE LASSO on ALL N observations using the
# cross-fitted (out-of-fold) nuisances in `nuis`. Combining out-of-fold
# stage-1 predictions with a FULL-SAMPLE stage-2 fit means the same
# observations inform both stages, which technically VIOLATES the theoretical
# sample-splitting / cross-fitting guarantees of Gao & Hastie (2022)
# Algorithm 1. It is kept deliberately as the documented contrast to the
# leakage-free learner_dina_crossfit() below.
# --------------------------------------------------------------------------
learner_dina_original <- function(X, Y, trt, nuis, X_eval = NULL,
                                  seed = NULL, nperm = 10) {
  X <- as.data.frame(X)
  for (cn in names(X)) if (is.character(X[[cn]])) X[[cn]] <- as.factor(X[[cn]])
  trt <- as.numeric(trt); Y <- as.numeric(Y)
  nu_hat <- nuis$e_hat * nuis$mu1_hat + (1 - nuis$e_hat) * nuis$mu0_hat
  Y_residual <- Y - nu_hat; W_centered <- trt - nuis$e_hat

  d <- .dina_design(X, W_centered)
  beta_final <- tryCatch({
    cv <- glmnet::cv.glmnet(d$WX_matrix, Y_residual, alpha = 1, intercept = FALSE,
                            penalty.factor = d$pen_factor, standardize = TRUE)
    as.numeric(coef(cv, s = "lambda.min"))[-1]
  }, error = function(e) {
    df <- as.data.frame(d$WX_matrix)
    f  <- as.formula(paste("yr ~ -1 +", paste(sprintf("`%s`", colnames(df)), collapse = " + ")))
    as.numeric(coef(lm(f, data = cbind(data.frame(yr = Y_residual), df))))
  })
  beta_final[is.na(beta_final)] <- 0
  names(beta_final) <- colnames(d$X_matrix)
  tau_hat  <- as.vector(d$X_matrix %*% beta_final)
  tau_eval <- if (!is.null(X_eval))
    .dina_eval_tau(X_eval, d$X_matrix, beta_final, mean(tau_hat, na.rm = TRUE)) else NULL

  hv <- .dina_het_vi(X, Y_residual, W_centered, tau_hat, nperm)
  # Top_Modifier (vi[1]) now uses the SAME definition as R-learner (LASSO):
  # |beta_j|*sd(x_j) from the stage-2 fit. cforest permimp kept as secondary.
  vi_coef <- .dina_coef_vi(beta_final, d$X_matrix)
  list(tau = tau_hat, tau_eval = tau_eval, method = "DINA Original",
       p_value = hv$het$p_max, p_value_max = hv$het$p_max, p_value_quad = hv$het$p_quad,
       vi = vi_coef, vi_coef = vi_coef, vi_cforest = hv$vi)
}

# --------------------------------------------------------------------------
# 4. DINA Cross-Fit — Gao & Hastie (2022) Algorithm 1, GENERALIZED to the
#    shared K=5 folds (port + generalization of dina_cross_fit..., :1050).
#
# Leakage-free DML1: iterate the EXACT folds in nuis$fold_ids. For each fold k,
# fit the stage-2 LASSO ENTIRELY on the rows in fold k — which use the
# out-of-fold nuisances (trained on folds != k), so stage-1 and stage-2 use
# disjoint data. Extract beta^(k); AVERAGE the K coefficient vectors.
# (Not forced to 2-fold and not collapsed — honors the K=5 same-footing nuis.)
# --------------------------------------------------------------------------
learner_dina_crossfit <- function(X, Y, trt, nuis, X_eval = NULL,
                                  seed = NULL, nperm = 10) {
  X <- as.data.frame(X)
  for (cn in names(X)) if (is.character(X[[cn]])) X[[cn]] <- as.factor(X[[cn]])
  trt <- as.numeric(trt); Y <- as.numeric(Y); n <- length(Y)
  nu_hat <- nuis$e_hat * nuis$mu1_hat + (1 - nuis$e_hat) * nuis$mu0_hat
  Y_residual <- Y - nu_hat; W_centered <- trt - nuis$e_hat

  d <- .dina_design(X, W_centered)
  coef_names <- colnames(d$X_matrix)
  fold_ids <- nuis$fold_ids; uf <- sort(unique(fold_ids))

  beta_list <- list()
  for (k_i in uf) {
    idx <- which(fold_ids == k_i)
    if (length(idx) < ncol(d$X_matrix) + 5) next   # too few rows to fit this fold
    bk <- tryCatch({
      cv <- glmnet::cv.glmnet(d$WX_matrix[idx, , drop = FALSE], Y_residual[idx],
                              alpha = 1, intercept = FALSE,
                              penalty.factor = d$pen_factor, standardize = TRUE)
      b <- as.numeric(coef(cv, s = "lambda.min"))[-1]; names(b) <- coef_names; b
    }, error = function(e) NULL)
    if (is.null(bk)) next
    bk[is.na(bk)] <- 0
    beta_list[[length(beta_list) + 1]] <- bk
  }

  if (length(beta_list) == 0) {
    return(list(tau = rep(NA_real_, n),
                tau_eval = if (!is.null(X_eval)) rep(NA_real_, nrow(as.data.frame(X_eval))) else NULL,
                method = "DINA Cross-Fit", p_value = NA, p_value_max = NA,
                p_value_quad = NA, vi = data.frame()))
  }
  # Average the per-fold beta vectors (DML1).
  beta_final <- Reduce("+", lapply(beta_list, function(b) b[coef_names])) / length(beta_list)
  names(beta_final) <- coef_names; beta_final[is.na(beta_final)] <- 0
  tau_hat  <- as.vector(d$X_matrix %*% beta_final)
  tau_eval <- if (!is.null(X_eval))
    .dina_eval_tau(X_eval, d$X_matrix, beta_final, mean(tau_hat, na.rm = TRUE)) else NULL

  hv <- .dina_het_vi(X, Y_residual, W_centered, tau_hat, nperm)
  # Same-footing Top_Modifier: coefficient VI (see DINA Original note).
  vi_coef <- .dina_coef_vi(beta_final, d$X_matrix)
  list(tau = tau_hat, tau_eval = tau_eval, method = "DINA Cross-Fit",
       p_value = hv$het$p_max, p_value_max = hv$het$p_max, p_value_quad = hv$het$p_quad,
       vi = vi_coef, vi_coef = vi_coef, vi_cforest = hv$vi,
       n_folds_used = length(beta_list))
}

# --------------------------------------------------------------------------
# 4b. DINA Gao (k2) — FAITHFUL Gao & Hastie (2022), SELF-CONTAINED.
#
# Two variants share an IDENTICAL self-contained stage 1 and differ ONLY in
# the stage-2 fit, so the GLM-vs-LASSO contrast is clean:
#   * DINA Gao (GLM k2)   — UNPENALIZED OLS stage 2 (the faithful default;
#                           paper §4.2 step (b) "MLE ... glm in R"). For
#                           Gaussian Y this IS the linear R-learner.
#   * DINA Gao (LASSO k2) — L1-penalized stage 2 (paper footnote: "penalties
#                           such as ridge and LASSO can be directly added").
#
# Neither consumes the shared 5-fold nuisances: Gao's Algorithm 1 COUPLES the
# nuisance split to the causal split (k=2; DINA.fit loop). For each half h we
# train stage-1 nuisances on the OTHER half, fit stage 2 on half h, and
# average the two coefficient vectors.
#
# Stage 1 (per half, Gao DINAEstimator.step1, Gaussian family):
#   nu(x) = E[Y|X] estimated DIRECTLY as the marginal mean ("we can directly
#   estimate nu(x) ... and avoid having to separately estimate eta0, eta1" —
#   paper §4.2), a(x) = e(x). Estimators mirror estimate_nuisances() exactly:
#   SuperLearner(SL.glmnet, SL.cforest)/method.CC_LS outcome, SL.glmnet/
#   method.NNLS2 propensity (estimated even in RCT, clipped [0.01,0.99]) —
#   SAME library, different split.
# --------------------------------------------------------------------------

# Shared self-contained stage 1 + assembly. `stage2` is a callback fitting one
# half's residualized design: (WX_sub, yr_sub, pen_factor, coef_names) -> beta.
.dina_gao_k2 <- function(X, Y, trt, X_eval, seed, nperm, method_label, stage2) {
  X <- as.data.frame(X)
  for (cn in names(X)) if (is.character(X[[cn]])) X[[cn]] <- as.factor(X[[cn]])
  trt <- as.numeric(trt); Y <- as.numeric(Y); n <- length(Y)
  Xn <- .coerce_numeric_df(X)   # same encoding as estimate_nuisances / R-learner

  outcome_lib <- list("SL.glmnet", "SL.cforest")
  prop_lib    <- "SL.glmnet"

  if (!is.null(seed)) set.seed(seed + 271828L)
  half <- sample(rep(1:2, length.out = n))

  # Per-row out-of-half nuisance predictions (row i gets the model trained on
  # the half it does NOT belong to) — used for stage 2 and the het test/VI.
  nu_hat <- rep(NA_real_, n); e_hat <- rep(NA_real_, n)
  for (h in 1:2) {
    test_idx  <- which(half == h)
    train_idx <- which(half != h)
    nu_hat[test_idx] <- tryCatch({
      f <- SuperLearner::SuperLearner(
        Y = Y[train_idx], X = Xn[train_idx, , drop = FALSE],
        SL.library = outcome_lib, family = "gaussian",
        method = "method.CC_LS", verbose = FALSE)
      as.numeric(predict(f, newdata = Xn[test_idx, , drop = FALSE])$pred)
    }, error = function(e) rep(mean(Y[train_idx]), length(test_idx)))
    e_hat[test_idx] <- tryCatch({
      f <- SuperLearner::SuperLearner(
        Y = trt[train_idx], X = Xn[train_idx, , drop = FALSE],
        SL.library = prop_lib, family = "binomial",
        method = "method.NNLS2", verbose = FALSE)
      as.numeric(predict(f, newdata = Xn[test_idx, , drop = FALSE])$pred)
    }, error = function(e) rep(mean(trt[train_idx]), length(test_idx)))
  }
  e_hat <- pmax(0.01, pmin(e_hat, 0.99))
  Y_residual <- Y - nu_hat; W_centered <- trt - e_hat

  d <- .dina_design(X, W_centered)
  coef_names <- colnames(d$X_matrix)

  beta_list <- list()
  for (h in 1:2) {
    idx <- which(half == h)
    bk <- tryCatch(
      stage2(d$WX_matrix[idx, , drop = FALSE], Y_residual[idx],
             d$pen_factor, coef_names),
      error = function(e) NULL)
    if (!is.null(bk)) { bk[is.na(bk)] <- 0; beta_list[[length(beta_list) + 1]] <- bk }
  }
  if (length(beta_list) == 0) {
    return(list(tau = rep(NA_real_, n),
                tau_eval = if (!is.null(X_eval)) rep(NA_real_, nrow(as.data.frame(X_eval))) else NULL,
                method = method_label, p_value = NA, p_value_max = NA,
                p_value_quad = NA, vi = data.frame()))
  }
  beta_final <- Reduce("+", lapply(beta_list, function(b) b[coef_names])) / length(beta_list)
  names(beta_final) <- coef_names; beta_final[is.na(beta_final)] <- 0
  tau_hat  <- as.vector(d$X_matrix %*% beta_final)
  tau_eval <- if (!is.null(X_eval))
    .dina_eval_tau(X_eval, d$X_matrix, beta_final, mean(tau_hat, na.rm = TRUE)) else NULL

  hv <- .dina_het_vi(X, Y_residual, W_centered, tau_hat, nperm)
  # Same-footing Top_Modifier: coefficient VI (see DINA Original note).
  vi_coef <- .dina_coef_vi(beta_final, d$X_matrix)
  list(tau = tau_hat, tau_eval = tau_eval, method = method_label,
       p_value = hv$het$p_max, p_value_max = hv$het$p_max, p_value_quad = hv$het$p_quad,
       vi = vi_coef, vi_coef = vi_coef, vi_cforest = hv$vi,
       n_folds_used = length(beta_list))
}

# Stage-2 callbacks ---------------------------------------------------------
# UNPENALIZED OLS (Gao DINAEstimator.step2 / paper §4.2 glm): faithful default.
.dina_gao_stage2_glm <- function(WX_sub, yr_sub, pen_factor, coef_names) {
  b <- lm.fit(WX_sub, yr_sub)$coefficients; names(b) <- coef_names; b
}
# L1-penalized (paper footnote): the same residualized design under cv.glmnet,
# matching the lasso config of DINA Original / Cross-Fit (lambda.min, the
# intercept column unpenalized via pen_factor). Falls back to OLS on error.
.dina_gao_stage2_lasso <- function(WX_sub, yr_sub, pen_factor, coef_names) {
  b <- tryCatch({
    cv <- glmnet::cv.glmnet(WX_sub, yr_sub, alpha = 1, intercept = FALSE,
                            penalty.factor = pen_factor, standardize = TRUE)
    as.numeric(coef(cv, s = "lambda.min"))[-1]
  }, error = function(e) lm.fit(WX_sub, yr_sub)$coefficients)
  names(b) <- coef_names; b
}

learner_dina_gao_glm <- function(X, Y, trt, nuis, X_eval = NULL,
                                 seed = NULL, nperm = 10)
  .dina_gao_k2(X, Y, trt, X_eval, seed, nperm,
               "DINA Gao (GLM k2)", .dina_gao_stage2_glm)

learner_dina_gao_lasso <- function(X, Y, trt, nuis, X_eval = NULL,
                                   seed = NULL, nperm = 10)
  .dina_gao_k2(X, Y, trt, X_eval, seed, nperm,
               "DINA Gao (LASSO k2)", .dina_gao_stage2_lasso)

# --------------------------------------------------------------------------
# 5. Dandl Extension — Dandl et al. (2024) model-based forest on residuals.
# Residualized base model + model4you::pmforest; OOB pmodel coefficients = tau.
# --------------------------------------------------------------------------
learner_dandl <- function(X, Y, trt, nuis, X_eval = NULL, seed = NULL, nperm = 10) {
  X <- as.data.frame(X); trt <- as.numeric(trt); Y <- as.numeric(Y)
  nu_hat <- nuis$e_hat * nuis$mu1_hat + (1 - nuis$e_hat) * nuis$mu0_hat
  Y_residual <- Y - nu_hat; W_centered <- trt - nuis$e_hat

  # Keep ALL covariates as partitioning variables. Model-based forests split on
  # factors natively (partykit), so categorical effect modifiers must NOT be
  # dropped — doing so blinds the forest to e.g. the benchtm continuous
  # scenarios whose true CATE depends on a categorical: scenario 3 pred =
  # (X14>0.25)&(X1=='N') and scenario 4 pred = (X14>0.3)|(X4=='Y'). Only
  # characters need coercion to factor; numeric columns stay numeric.
  X_part <- X
  for (cn in names(X_part)) if (is.character(X_part[[cn]])) X_part[[cn]] <- as.factor(X_part[[cn]])
  xnames   <- names(X_part)
  est_data <- data.frame(Y_residual = Y_residual, W_centered = W_centered,
                         X_part, check.names = FALSE)
  # Partition only on the covariates (exclude the base-model variables), so
  # mtry = P refers to the covariate dimension exactly.
  zform    <- stats::as.formula(paste("~", paste(sprintf("`%s`", xnames), collapse = " + ")))

  base_model   <- lm(Y_residual ~ W_centered, data = est_data)
  # Dandl et al. (2024) / htesim Section 7 hyperparameters for the continuous
  # case: M = 500 trees, minimum node size = 14, mtry = P (every covariate is a
  # split candidate — the defining feature of model-based forests), subsampling.
  ctrl         <- partykit::ctree_control(teststat = "quad", testtype = "Univ",
                                          mincriterion = 0, minbucket = 14L,
                                          saveinfo = FALSE, lookahead = TRUE)
  forest_model <- model4you::pmforest(base_model, data = est_data, zformula = zform,
                                      ntree = 500L, mtry = length(xnames),
                                      perturb = list(replace = FALSE, fraction = 0.632),
                                      control = ctrl)
  tau_hat      <- model4you::pmodel(forest_model, OOB = TRUE)[, "W_centered"]

  tau_eval <- NULL
  if (!is.null(X_eval)) {
    Xe <- as.data.frame(X_eval)
    for (cn in names(Xe)) if (is.character(Xe[[cn]])) Xe[[cn]] <- as.factor(Xe[[cn]])
    # Align factor levels to the training data so new obs can be placed in leaves.
    for (cn in xnames) if (is.factor(est_data[[cn]]) && cn %in% names(Xe))
      Xe[[cn]] <- factor(as.character(Xe[[cn]]), levels = levels(est_data[[cn]]))
    nd <- data.frame(Y_residual = 0, W_centered = 0,
                     Xe[, xnames, drop = FALSE], check.names = FALSE)
    tau_eval <- tryCatch(
      as.vector(model4you::pmodel(forest_model, newdata = nd)[, "W_centered"]),
      error = function(e) rep(mean(tau_hat, na.rm = TRUE), nrow(Xe)))
  }

  hv <- .dina_het_vi(X_part, Y_residual, W_centered, tau_hat, nperm)
  list(tau = tau_hat, tau_eval = tau_eval, method = "Dandl Extension",
       p_value = hv$het$p_max, p_value_max = hv$het$p_max, p_value_quad = hv$het$p_quad,
       vi = hv$vi)
}

# --------------------------------------------------------------------------
# 6. R-learner (LASSO) — Nie & Wager (2021) via rlearner::rlasso.
# rlearner::rlasso fed the SHARED nuisances (p_hat = e_hat, m_hat = composite),
# so it isolates the stage-2 estimator. Global test on the R-learner pseudo-tau.
# --------------------------------------------------------------------------
learner_rlearner_lasso <- function(X, Y, trt, nuis, X_eval = NULL,
                                   seed = NULL, nperm = 10) {
  if (!requireNamespace("rlearner", quietly = TRUE)) return(NULL)
  Xc <- .coerce_numeric_df(X); X_matrix <- as.matrix(Xc)
  m_hat <- nuis$e_hat * nuis$mu1_hat + (1 - nuis$e_hat) * nuis$mu0_hat
  rl_fit <- rlearner::rlasso(X_matrix, as.numeric(trt), as.numeric(Y),
                             p_hat = nuis$e_hat, m_hat = m_hat)
  tau_hat <- as.vector(predict(rl_fit, X_matrix))

  tau_eval <- NULL
  if (!is.null(X_eval)) {
    Xe_m <- as.matrix(.coerce_numeric_df(X_eval))
    tau_eval <- tryCatch(as.vector(predict(rl_fit, Xe_m)),
                         error = function(e) rep(mean(tau_hat, na.rm = TRUE), nrow(Xe_m)))
  }

  het <- global_heterogeneity_test(
    .r_learner_pseudo(Y, trt, nuis$e_hat, nuis$mu0_hat, nuis$mu1_hat), X)
  vi_coef <- tryCatch({
    beta_vec <- as.numeric(rl_fit$tau_beta)
    if (length(beta_vec) != ncol(X_matrix) + 1L) stop("tau_beta length mismatch")
    beta_vec <- beta_vec[-1]; names(beta_vec) <- colnames(X_matrix)
    # rlasso standardizes X internally (caret center+scale), so tau_beta is
    # already per-SD; .dina_coef_vi's |beta_j| * sd(x_j) would weight by ~sd^2.
    imp <- abs(beta_vec); imp[!is.finite(imp)] <- 0
    vi <- data.frame(variable = colnames(X_matrix), importance = as.numeric(imp),
                     stringsAsFactors = FALSE)
    vi[order(-vi$importance, sample(nrow(vi))), ]  # random tie-break
  }, error = function(e) NULL)

  list(tau = tau_hat, tau_eval = tau_eval, method = "R-learner (LASSO)",
       p_value = het$p_max, p_value_max = het$p_max, p_value_quad = het$p_quad,
       vi = vi_coef)
}

# --------------------------------------------------------------------------
# 7. TMLE (Li 2026 VTE) — CATE-targeted TMLE for the variance of the ITE.
# CATE-targeted TMLE seeded from the SHARED nuisances; returns the VTE point
# estimate + EIC-based CI (heterogeneity "detected" if CI excludes 0).
# --------------------------------------------------------------------------
learner_tmle_vte <- function(X, Y, trt, nuis, X_eval = NULL,
                             seed = NULL, nperm = 10, max_iter = 25) {
  trt <- as.numeric(trt); Y <- as.numeric(Y); n <- length(Y)
  e_hat <- pmax(0.01, pmin(nuis$e_hat, 0.99))
  mu0_curr <- as.numeric(nuis$mu0_hat); mu1_curr <- as.numeric(nuis$mu1_hat)
  tau_curr <- mu1_curr - mu0_curr; Q_bar <- ifelse(trt == 1, mu1_curr, mu0_curr)

  converged <- FALSE; iters <- 0L
  for (it in seq_len(max_iter)) {
    iters <- it; tau_mean <- mean(tau_curr)
    H <- 2 * (tau_curr - tau_mean) * (2 * trt - 1) / ifelse(trt == 1, e_hat, 1 - e_hat)
    D <- H * (Y - Q_bar); sigma_D <- sd(D)
    tol <- if (is.finite(sigma_D) && sigma_D > 0) sigma_D / (sqrt(n) * log(n)) else 1e-4
    if (abs(mean(D)) <= tol) { converged <- TRUE; break }
    eps_fit <- tryCatch(glm(Y ~ -1 + H + offset(Q_bar), family = gaussian()),
                        error = function(e) NULL)
    eps_it <- if (is.null(eps_fit)) 0 else { v <- unname(coef(eps_fit)["H"])
                                             if (is.na(v) || !is.finite(v)) 0 else v }
    if (abs(eps_it) > 0.5) eps_it <- sign(eps_it) * 0.5
    H_treated <-  2 * (tau_curr - tau_mean) / e_hat
    H_control <- -2 * (tau_curr - tau_mean) / (1 - e_hat)
    mu1_curr <- mu1_curr + eps_it * H_treated; mu0_curr <- mu0_curr + eps_it * H_control
    tau_curr <- mu1_curr - mu0_curr; Q_bar <- ifelse(trt == 1, mu1_curr, mu0_curr)
  }

  tau_tmle <- as.vector(tau_curr); tau_mean_fin <- mean(tau_tmle)
  vte <- mean((tau_tmle - tau_mean_fin)^2)
  H_fin <- 2 * (tau_tmle - tau_mean_fin) * (2 * trt - 1) / ifelse(trt == 1, e_hat, 1 - e_hat)
  IC <- H_fin * (Y - Q_bar) + (tau_tmle - tau_mean_fin)^2 - vte
  vte_se <- sqrt(var(IC) / n)
  vte_ci_lo <- vte - 1.96 * vte_se; vte_ci_hi <- vte + 1.96 * vte_se

  pseudo_dr <- tau_tmle + (trt - e_hat) / (e_hat * (1 - e_hat)) * (Y - Q_bar)
  het <- global_heterogeneity_test(pseudo_dr, X)

  list(tau = tau_tmle, tau_eval = if (!is.null(X_eval)) rep(NA_real_, nrow(as.data.frame(X_eval))) else NULL,
       method = "TMLE (Li 2026 VTE)",
       vte = vte, vte_se = vte_se, vte_ci_lo = vte_ci_lo, vte_ci_hi = vte_ci_hi,
       vte_detected = if (is.finite(vte_ci_lo)) vte_ci_lo > 0 else NA,
       n_iter = iters, converged = converged,
       p_value = het$p_max, p_value_max = het$p_max, p_value_quad = het$p_quad)
}

# --------------------------------------------------------------------------
# 8/9. WATCH score-residual (Residual LASSO / Residual Risk).
# After the WATCH score-residual approach (exploreTEH, Sechidis et al.).
#
# Model-based heterogeneity TEST + VI. By design it builds its OWN prognostic
# risk score (riskPred via cv.glmnet; alpha = 1 LASSO / alpha = 0 ridge) — it
# does NOT use the outcome nuisances. It consumes nuis$e_hat only for the
# centered-treatment convention (unused for the score itself; exploreTEH_mb
# centers on mean(trt)). A parametric CATE proxy keeps tau / Somers_D populated.
# --------------------------------------------------------------------------
.watch_lasso_select <- function(X_mat, Y, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  cvfit <- tryCatch(glmnet::cv.glmnet(X_mat, Y, alpha = 1, nfolds = 10, family = "gaussian"),
                    error = function(e) NULL)
  if (is.null(cvfit)) return(seq_len(min(5L, ncol(X_mat))))
  coefs <- as.vector(coef(cvfit, s = "lambda.min"))[-1]
  sel <- which(abs(coefs) > 1e-8)
  if (length(sel) == 0) sel <- order(abs(coefs), decreasing = TRUE)[1]
  sel
}
.watch_risk_pred <- function(X_mat, Y, X_eval_mat = NULL, alpha = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  cvfit <- tryCatch(glmnet::cv.glmnet(X_mat, Y, alpha = alpha, nfolds = 10, family = "gaussian"),
                    error = function(e) NULL)
  if (is.null(cvfit)) return(list(
    train = rep(mean(Y, na.rm = TRUE), nrow(X_mat)),
    eval  = if (!is.null(X_eval_mat)) rep(mean(Y, na.rm = TRUE), nrow(X_eval_mat)) else NULL))
  list(train = as.vector(predict(cvfit, newx = X_mat, s = "lambda.min")),
       eval  = if (!is.null(X_eval_mat)) as.vector(predict(cvfit, newx = X_eval_mat, s = "lambda.min")) else NULL)
}
.watch_cate_proxy <- function(Y, z_tilde, prog_train, X_pred_train,
                              prog_eval = NULL, X_pred_eval = NULL) {
  Y <- as.numeric(Y); z_tilde <- as.numeric(z_tilde)
  q_prog <- ncol(prog_train); q_pred <- ncol(X_pred_train)
  inter_train <- if (q_pred > 0) X_pred_train * z_tilde else matrix(0, nrow = length(Y), ncol = 0)
  D_train <- cbind(rep(1, length(Y)),
                   if (q_prog > 0) prog_train else matrix(0, nrow = length(Y), ncol = 0),
                   z_tilde, inter_train)
  colnames(D_train) <- c("(Intercept)", if (q_prog > 0) paste0("prog_", seq_len(q_prog)),
                         "z", if (q_pred > 0) paste0("zX_", seq_len(q_pred)))
  fit <- tryCatch(qr.solve(D_train, Y), error = function(e) NULL)
  if (is.null(fit)) return(list(tau = rep(0, length(Y)), tau_eval = NULL))
  names(fit) <- colnames(D_train)
  delta_hat <- as.numeric(fit["z"])
  gamma_hat <- if (q_pred > 0) as.numeric(fit[paste0("zX_", seq_len(q_pred))]) else numeric(0)
  gamma_hat[!is.finite(gamma_hat)] <- 0
  predict_tau <- function(B) if (q_pred == 0) rep(delta_hat, nrow(B)) else
    delta_hat + as.vector(B %*% gamma_hat)
  list(tau = predict_tau(X_pred_train),
       tau_eval = if (!is.null(X_pred_eval)) predict_tau(X_pred_eval) else NULL)
}
.watch_score_residual_core <- function(X, Y, trt, alpha, method_label,
                                       nperm, X_eval, seed, e_hat = NULL) {
  X <- .coerce_numeric_df(X); Y_n <- as.numeric(Y); trt_n <- as.numeric(trt)
  X_mat <- as.matrix(X)
  Xe_df  <- if (!is.null(X_eval)) .coerce_numeric_df(X_eval) else NULL
  Xe_mat <- if (!is.null(Xe_df)) as.matrix(Xe_df) else NULL

  risk <- .watch_risk_pred(X_mat, Y_n, Xe_mat, alpha = alpha, seed = seed)
  trt_centered <- trt_n - mean(trt_n)
  fit_homo <- tryCatch(stats::glm(Y_n ~ trt_centered + risk$train, family = stats::gaussian()),
                       error = function(e) NULL)
  if (is.null(fit_homo)) return(NULL)
  s_i <- tryCatch(as.numeric(sandwich::estfun(fit_homo)[, "trt_centered"]),
                  error = function(e) NULL)
  if (is.null(s_i)) return(NULL)

  het <- global_heterogeneity_test(s_i, X)
  vi <- tryCatch({
    cf <- party::cforest(s_i ~ ., data = data.frame(X, s_i = s_i),
                         control = party::cforest_unbiased(mtry = 5, ntree = 500))
    v  <- permimp::permimp(cf, conditional = FALSE, nperm = nperm)
    df <- data.frame(variable = names(v$values), importance = v$values, stringsAsFactors = FALSE)
    df[order(-df$importance), ]
  }, error = function(e)
    data.frame(variable = names(X), importance = rep(0, ncol(X)), stringsAsFactors = FALSE))

  sel <- .watch_lasso_select(X_mat, Y_n, seed = seed)
  cate <- .watch_cate_proxy(
    Y = Y_n, z_tilde = trt_centered,
    prog_train = matrix(risk$train, ncol = 1, dimnames = list(NULL, "risk")),
    X_pred_train = X_mat[, sel, drop = FALSE],
    prog_eval = if (!is.null(risk$eval)) matrix(risk$eval, ncol = 1) else NULL,
    X_pred_eval = if (!is.null(Xe_mat)) Xe_mat[, sel, drop = FALSE] else NULL)

  list(tau = cate$tau, tau_eval = cate$tau_eval, method = method_label,
       p_value = het$p_max, p_value_max = het$p_max, p_value_quad = het$p_quad, vi = vi)
}
learner_watch_lasso <- function(X, Y, trt, nuis, X_eval = NULL, seed = NULL, nperm = 10) {
  .watch_score_residual_core(X, Y, trt, alpha = 1, method_label = "WATCH (Residual LASSO)",
                             nperm = nperm, X_eval = X_eval, seed = seed, e_hat = nuis$e_hat)
}
learner_watch_risk <- function(X, Y, trt, nuis, X_eval = NULL, seed = NULL, nperm = 10) {
  .watch_score_residual_core(X, Y, trt, alpha = 0, method_label = "WATCH (Residual Risk)",
                             nperm = nperm, X_eval = X_eval, seed = seed, e_hat = nuis$e_hat)
}

# ==========================================================================
# §3z. DINA CENTERED — rlasso-parity stage-2 design. Shares the SAME
# nuisances as DINA Original; differs ONLY in the stage-2 construction:
#
#   * X is center+scaled BEFORE forming the (W - e_hat) * [1, X] interaction
#     (rlasso uses caret center+scale first), so cv.glmnet runs with
#     standardize=FALSE and a FREE (unpenalized) intercept.
#   * Pooled stage 2 on all n (DML2 / rlasso architecture); the nuisances are
#     already out-of-fold for every row, so the quasi-oracle guarantee holds.
#
# Why: with the RAW-X design of DINA Original the penalized resW.Xj columns are
# correlated with the unpenalized resW (ATE) column, leaking ATE signal into
# modifier selection — the cause of Original's weak S3/S4 top-modifier hits.
# Centering removes that coupling. In our tuning experiments DINA Centered
# matched the R-learner on all three WATCH objectives and significantly beat
# DINA Original (Obj-2 beta=2 hit 0.635 vs 0.515; MSE -0.009, p=3e-11).
#
# DINA Centered Expanded additionally appends quartile-threshold indicators
# I(x > q25/q50/q75) for continuous covariates (targets the benchtm S3/S4
# subgroup pred functions) and aggregates Top_Modifier back to parent
# covariates — best CATE accuracy of the tested variants.
#
# Both keep DINA Original / DINA Cross-Fit untouched as faithful-replica
# controls.
# ==========================================================================

# Centered design from a numeric feature matrix.
.dina_design_centered_mat <- function(feats, W_centered) {
  feats <- as.matrix(feats)
  mu  <- colMeans(feats)
  sdv <- apply(feats, 2, stats::sd)
  sdv[!is.finite(sdv) | sdv < 1e-12] <- 1
  Xs  <- sweep(sweep(feats, 2, mu, "-"), 2, sdv, "/")
  X_matrix <- cbind(`(Intercept)` = 1, Xs)
  WX <- W_centered * X_matrix
  colnames(WX) <- c("resW", paste0("resW.", colnames(Xs)))
  list(X_matrix = X_matrix, WX_matrix = WX,
       pen_factor = c(0, rep(1, ncol(Xs))),
       center = mu, scale = sdv, feat_names = colnames(Xs))
}

# Threshold-feature expansion for continuous covariates (quartile indicators).
# `thresholds` carries the TRAINING quantiles to the eval set.
.dina_expand_features <- function(X, thresholds = NULL) {
  Xdf <- as.data.frame(X)
  Xn  <- .coerce_numeric_df(Xdf)
  if (is.null(thresholds)) {
    cont <- names(Xdf)[vapply(Xdf, function(cl)
      is.numeric(cl) && length(unique(cl)) > 5, logical(1))]
    thresholds <- lapply(Xn[cont], function(v)
      unique(stats::quantile(v, c(.25, .5, .75), na.rm = TRUE, names = FALSE)))
    names(thresholds) <- cont
  }
  feats <- Xn
  for (v in names(thresholds)) {
    if (!v %in% names(Xn)) next
    qs <- thresholds[[v]]
    for (qi in seq_along(qs))
      feats[[paste0(v, "..gt", qi)]] <- as.numeric(Xn[[v]] > qs[qi])
  }
  list(feats = as.matrix(feats), thresholds = thresholds)
}

# Parent-covariate VI for expanded designs: sum |beta_f|*sd(f) over the
# features derived from each covariate, then rank parents (random tie-break).
.dina_group_coef_vi <- function(beta_final, X_matrix) {
  vi_f <- .dina_coef_vi(beta_final, X_matrix)
  if (nrow(vi_f) == 0) return(vi_f)
  vi_f$parent <- sub("\\.\\.gt[0-9]+$", "", vi_f$variable)
  agg <- stats::aggregate(importance ~ parent, data = vi_f, FUN = sum)
  names(agg) <- c("variable", "importance")
  agg[order(-agg$importance, sample(nrow(agg))), ]
}

# Het test only (no cforest permimp) — same pseudo-outcome as .dina_het_vi.
.dina_het_only <- function(X, Y_residual, W_centered) {
  pseudo_tau <- rep(NA_real_, length(W_centered)); nz <- abs(W_centered) > 0.01
  pseudo_tau[nz] <- Y_residual[nz] / W_centered[nz]
  ok <- !is.na(pseudo_tau) & is.finite(pseudo_tau)
  if (sum(ok) > 10) global_heterogeneity_test(pseudo_tau[ok], X[ok, , drop = FALSE])
  else list(p_max = NA_real_, p_quad = NA_real_)
}

# Centered, pooled stage-2 core (design = "centered" | "expanded").
.dina_centered_core <- function(X, Y, trt, nuis, X_eval, seed, nperm,
                                method_label, design = "centered") {
  X <- as.data.frame(X)
  for (cn in names(X)) if (is.character(X[[cn]])) X[[cn]] <- as.factor(X[[cn]])
  trt <- as.numeric(trt); Y <- as.numeric(Y); n <- length(Y)
  nu_hat <- nuis$e_hat * nuis$mu1_hat + (1 - nuis$e_hat) * nuis$mu0_hat
  Y_residual <- Y - nu_hat; W_centered <- trt - nuis$e_hat

  expanded <- NULL
  feats <- if (design == "expanded") {
    expanded <- .dina_expand_features(X); expanded$feats
  } else as.matrix(.coerce_numeric_df(X))
  d <- .dina_design_centered_mat(feats, W_centered)
  coef_names <- colnames(d$X_matrix)

  beta_final <- tryCatch({
    cv <- glmnet::cv.glmnet(d$WX_matrix, Y_residual, alpha = 1, intercept = TRUE,
                            penalty.factor = d$pen_factor, standardize = FALSE)
    b <- as.numeric(coef(cv, s = "lambda.min"))[-1]; names(b) <- coef_names; b
  }, error = function(e) { b <- rep(0, length(coef_names)); names(b) <- coef_names; b })
  beta_final[is.na(beta_final)] <- 0
  tau_hat <- as.vector(d$X_matrix %*% beta_final)

  tau_eval <- NULL
  if (!is.null(X_eval)) {
    tau_eval <- tryCatch({
      Xe <- as.data.frame(X_eval)
      feats_e <- if (design == "expanded")
        .dina_expand_features(Xe, thresholds = expanded$thresholds)$feats
      else as.matrix(.coerce_numeric_df(Xe))
      feats_e <- feats_e[, d$feat_names, drop = FALSE]
      Xs_e <- sweep(sweep(feats_e, 2, d$center, "-"), 2, d$scale, "/")
      as.vector(cbind(1, Xs_e) %*% beta_final)
    }, error = function(e) rep(mean(tau_hat, na.rm = TRUE), nrow(as.data.frame(X_eval))))
  }

  het <- .dina_het_only(X, Y_residual, W_centered)
  vi_coef <- if (design == "expanded") .dina_group_coef_vi(beta_final, d$X_matrix)
             else .dina_coef_vi(beta_final, d$X_matrix)
  list(tau = tau_hat, tau_eval = tau_eval, method = method_label,
       p_value = het$p_max, p_value_max = het$p_max, p_value_quad = het$p_quad,
       vi = vi_coef, vi_coef = vi_coef)
}

learner_dina_centered <- function(X, Y, trt, nuis, X_eval = NULL, seed = NULL, nperm = 10)
  .dina_centered_core(X, Y, trt, nuis, X_eval, seed, nperm, "DINA Centered",
                      design = "centered")
learner_dina_centered_exp <- function(X, Y, trt, nuis, X_eval = NULL, seed = NULL, nperm = 10)
  .dina_centered_core(X, Y, trt, nuis, X_eval, seed, nperm, "DINA Centered Expanded",
                      design = "expanded")


# ==========================================================================
# §4. DRIVER:  compare_meta_learners()
#
# Estimates the shared nuisances ONCE, dispatches every requested learner on
# those SAME nuisances, and (if an eval set is given) computes the Gao-2024
# pairwise relative-error EIF between methods.
# ==========================================================================

# Registry: method label -> learner function. The reference is listed first.
METHOD_REGISTRY <- list(
  "DR-learner (2-Model 2SL)" = learner_dr_2model_2sl,   # REFERENCE
  "DR-learner (Kennedy)"     = learner_dr_kennedy,
  "DINA Original"            = learner_dina_original,
  "DINA Cross-Fit"           = learner_dina_crossfit,
  "DINA Gao (GLM k2)"        = learner_dina_gao_glm,
  "DINA Gao (LASSO k2)"      = learner_dina_gao_lasso,
  "Dandl Extension"          = learner_dandl,
  "R-learner (LASSO)"        = learner_rlearner_lasso,
  "TMLE (Li 2026 VTE)"       = learner_tmle_vte,
  "WATCH (Residual LASSO)"   = learner_watch_lasso,
  "WATCH (Residual Risk)"    = learner_watch_risk,
  "DINA Centered"            = learner_dina_centered,            # §3z (promoted 2026-06-12)
  "DINA Centered Expanded"   = learner_dina_centered_exp         # §3z (promoted 2026-06-12)
)

compare_meta_learners <- function(X, Y, trt,
                                  methods = names(METHOD_REGISTRY),
                                  tau_true = NULL,
                                  X_eval = NULL, Y_eval = NULL, trt_eval = NULL,
                                  seed = NULL, treatment_scenario = "RCT",
                                  light_config = FALSE, nperm = 10, verbose = TRUE,
                                  compute_vim = FALSE, vim_subsets = NULL) {
  X <- as.data.frame(X)
  unknown <- setdiff(methods, names(METHOD_REGISTRY))
  if (length(unknown) > 0) stop("Unknown method(s): ", paste(unknown, collapse = ", "))

  # --- Stage 1: shared nuisances (ONCE) ----------------------------------
  if (verbose) cat("== Estimating shared nuisances (Kostas 2-Model 2SL config) ==\n")
  nuis <- estimate_nuisances(X, Y, trt, seed = seed,
                             treatment_scenario = treatment_scenario,
                             light_config = light_config, verbose = verbose)

  # --- Eval-set nuisances for the Gao EIF (out-of-sample) ----------------
  nuis_eval <- NULL
  have_eval <- !is.null(X_eval) && !is.null(Y_eval) && !is.null(trt_eval)
  if (have_eval) {
    nuis_eval <- tryCatch(
      estimate_nuisances_eval(X, Y, trt, X_eval, seed = seed,
                              treatment_scenario = treatment_scenario,
                              light_config = light_config),
      error = function(e) { if (verbose) cat("  eval nuisance fit failed:", e$message, "\n"); NULL })
    if (is.null(nuis_eval)) have_eval <- FALSE
  }
  X_eval_dispatch <- if (!is.null(X_eval)) X_eval else NULL

  # --- Stage 2: dispatch each learner on the SAME nuisances --------------
  results <- list(); timing <- list()
  for (mname in methods) {
    fn <- METHOD_REGISTRY[[mname]]
    if (verbose) cat(sprintf("  -> %s\n", mname))
    t0 <- proc.time()[["elapsed"]]
    out <- tryCatch(fn(X, Y, trt, nuis, X_eval = X_eval_dispatch, seed = seed, nperm = nperm),
                    error = function(e) { if (verbose) cat("     failed:", e$message, "\n"); NULL })
    timing[[mname]] <- proc.time()[["elapsed"]] - t0
    if (!is.null(out)) results[[mname]] <- out
  }

  # --- Metrics table -----------------------------------------------------
  metrics_rows <- list()
  for (mname in names(results)) {
    r <- results[[mname]]
    em <- if (!is.null(tau_true)) evaluate_method(r$tau, tau_true) else
      list(MSE = NA, Bias = NA, Variance = NA, Correlation = NA, MAE = NA, RMSE = NA, Somers_D = NA)
    p_max  <- .num_or_na(if (!is.null(r$p_value_max)) r$p_value_max else r$p_value)
    p_quad <- .num_or_na(r$p_value_quad)
    metrics_rows[[mname]] <- data.frame(
      Method = r$method, Treatment_Scenario = treatment_scenario,
      MSE = em$MSE, Bias = em$Bias, Variance = em$Variance, Correlation = em$Correlation,
      MAE = em$MAE, RMSE = em$RMSE, Somers_D = em$Somers_D,
      Het_P_Value_Max = p_max, Het_P_Value_Quad = p_quad,
      Het_Detected_Max  = if (is.na(p_max))  NA else p_max  < 0.05,
      Het_Detected_Quad = if (is.na(p_quad)) NA else p_quad < 0.05,
      Top_Modifier = if (!is.null(r$vi) && nrow(r$vi) > 0) r$vi$variable[1] else NA,
      Time_Sec = .num_or_na(timing[[mname]]),
      VTE_Estimate = .num_or_na(r$vte),  VTE_SE = .num_or_na(r$vte_se),
      VTE_CI_Lo = .num_or_na(r$vte_ci_lo), VTE_CI_Hi = .num_or_na(r$vte_ci_hi),
      VTE_Detected = if (!is.null(r$vte_detected)) r$vte_detected else NA,
      stringsAsFactors = FALSE)
  }
  metrics <- if (length(metrics_rows) > 0) do.call(rbind, metrics_rows) else NULL
  if (!is.null(metrics)) rownames(metrics) <- NULL

  # --- Pairwise relative-error EIF (Gao 2024) ----------------------------
  pairwise <- NULL
  if (have_eval) {
    valid_m <- names(results)[vapply(names(results), function(m) {
      te <- results[[m]]$tau_eval
      !is.null(te) && is.numeric(te) && length(te) == length(Y_eval) && any(is.finite(te))
    }, logical(1))]
    rows <- list()
    if (length(valid_m) >= 2) {
      for (i in seq_along(valid_m)) for (j in seq_along(valid_m)) {
        if (j <= i) next
        m1 <- valid_m[i]; m2 <- valid_m[j]
        reif <- tryCatch(compute_relative_error_eif(
          tau1 = results[[m1]]$tau_eval, tau2 = results[[m2]]$tau_eval,
          Y = Y_eval, W = trt_eval,
          e_hat = nuis_eval$e_hat, mu0_hat = nuis_eval$mu0_hat, mu1_hat = nuis_eval$mu1_hat),
          error = function(e) list(delta = NA_real_, se = NA_real_, ci_lo = NA_real_,
                                   ci_hi = NA_real_, ci_excludes_zero = NA, n = 0L))
        winner <- if (is.na(reif$ci_excludes_zero) || !reif$ci_excludes_zero) "Tie" else
          if (reif$delta < 0) m1 else m2
        rows[[length(rows) + 1]] <- data.frame(
          Method1 = m1, Method2 = m2, Delta = reif$delta, SE = reif$se,
          CI_Lo = reif$ci_lo, CI_Hi = reif$ci_hi,
          CI_Excludes_Zero = reif$ci_excludes_zero, Winner = winner, N_Obs = reif$n,
          stringsAsFactors = FALSE)
      }
    }
    if (length(rows) > 0) { pairwise <- do.call(rbind, rows); rownames(pairwise) <- NULL }
  }

  # --- TE-VIM TMLE (Li 2026 Psi_2 / Psi_3) on the SAME shared nuisances ----
  # Per-covariate leave-one-out VIM; tau^0 is the S-learner mu1_hat - mu0_hat
  # derived inside compute_te_vim() (no dependency on any CATE learner above).
  tevims <- NULL
  if (compute_vim && exists("compute_te_vim", mode = "function")) {
    if (verbose) cat("  -> TE-VIM TMLE (Li 2026 Psi_2/Psi_3, per-covariate LOO)\n")
    tevims <- tryCatch(
      compute_te_vim(X, Y, trt, nuis, subsets = vim_subsets, verbose = FALSE),
      error = function(e) { if (verbose) cat("     TE-VIM failed:", e$message, "\n"); NULL })
  }

  list(results = results, metrics = metrics, pairwise_relative_errors = pairwise,
       tevims_scores = tevims, nuisance = nuis)
}


# ==========================================================================
# §5. GUARDED DEMO / SMOKE TEST
#
# Runs only when the file is executed directly from the repo root
# (Rscript R/meta_learners.R) and COMP_CLEAN_NO_DEMO is not set. Sources
# R/data_generation.R to reuse generate_simulation_data().
# ==========================================================================
if (sys.nframe() == 0 && !exists("COMP_CLEAN_NO_DEMO")) {
  message("== meta_learners demo ==")
  if (file.exists("R/data_generation.R")) {
    source("R/data_generation.R")

    seed_tr <- 12345L; seed_ev <- seed_tr + 999999L
    n_demo <- 300L; p_demo <- 30L; scen <- 1L; beta_lvl <- 1L

    dtrain <- generate_simulation_data(n = n_demo, p = p_demo, scenario = scen,
                                       beta_level = beta_lvl, seed = seed_tr,
                                       treatment_scenario = "RCT")
    deval  <- generate_simulation_data(n = n_demo, p = p_demo, scenario = scen,
                                       beta_level = beta_lvl, seed = seed_ev,
                                       treatment_scenario = "RCT")

    X <- dtrain %>% dplyr::select(dplyr::starts_with("X")); Y <- dtrain$Y; trt <- dtrain$trt
    tau_true <- dtrain$trt_effect
    X_eval <- deval %>% dplyr::select(dplyr::starts_with("X"))
    Y_eval <- deval$Y; trt_eval <- deval$trt

    res <- compare_meta_learners(
      X, Y, trt, tau_true = tau_true,
      X_eval = X_eval, Y_eval = Y_eval, trt_eval = trt_eval,
      seed = seed_tr, treatment_scenario = "RCT",
      light_config = TRUE, nperm = 5, verbose = TRUE,
      compute_vim = TRUE)

    cat("\n===== METRICS =====\n"); print(res$metrics)
    cat("\n===== PAIRWISE RELATIVE-ERROR (Gao 2024 EIF) =====\n")
    print(res$pairwise_relative_errors)
    cat("\n===== TE-VIM TMLE (Li 2026 Psi_2/Psi_3) — top 8 by VIMa =====\n")
    if (!is.null(res$tevims_scores)) {
      tv <- res$tevims_scores[order(-res$tevims_scores$vima), ]
      print(utils::head(tv[, c("variable","vima","vima_ci_l","vima_ci_u","vimb","vte")], 8))
      cat(sprintf("True modifier(s) for scenario %d should rank at/near the top.\n", scen))
    } else cat("(no TE-VIM scores returned)\n")
  } else {
    message("R/data_generation.R not found — run the demo from the repo root; ",
            "or call compare_meta_learners(X, Y, trt, ...) directly.")
  }
}
