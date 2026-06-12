############################################################################
# R/te_vim.R
#
# Faithful port of the TE-VIM TMLE from HaodongL/te_vim (the code companion to
# Li, Hubbard, Hines, Storås, Kvist & van der Laan, "Targeted Learning on a
# Variable Importance Measure for Heterogeneous Treatment Effects", 2026).
#
# WHAT THIS IMPLEMENTS
# --------------------
# The paper's THREE target statistical parameters (section 2.4):
#   Psi_1 = var(tau(W))                       VTE  (companion, Levy 2021)
#   Psi_2 = var(tau(W)) - var(tau_s(W))       VIMa (the paper's contribution)
#   Psi_3 = Psi_2 / Psi_1                      VIMb (scaled VIM)
# where tau_s(W) = E(tau(W) | W_{-s}) is the CATE projected onto the covariates
# EXCLUDING the subset s. Psi_2 measures the heterogeneity explained by s.
#
# SOURCE MAP (HaodongL/te_vim @ master)
#   R/est_function/vim.R       -> TMLE_VIM (Psi_2), TMLE_VTE (Psi_1), bound
#   R/est_function/fit_para.R  -> fit_para / fit_x: build tau, tau_s, gamma_s
# The two iterative TMLE targeting loops below (TMLE_VIM / TMLE_VTE) are ported
# VERBATIM from vim.R: same variable names, the same sign(PnD) coordinate-step
# update rule, the same sigma/(sqrt(N) log N) stopping criterion, the same
# glm(tau*^2 ~ 1, offset = gamma_s_0) update for gamma_s, and the influence
# curve evaluated at the INITIAL fits.
#
# FOUR REVIEWED ADAPTATIONS for this harness (flagged inline as ## [A1]-[A4]):
#   [A1] init: tau^0 is the S-LEARNER mu1_hat - mu0_hat from the shared
#        nuisances (NOT a DR-learner). This keeps tau^(i) = Q1^(i) - Q0^(i)
#        consistent at EVERY iteration including i = 0, matching the paper's
#        Algorithm 1 (where tau is derived from Q-bar). The repo's dr = TRUE
#        default seeds tau^0 with a DR fit that iteration 1 immediately
#        overwrites with the fluctuated S-learner tau -> an init inconsistency
#        we deliberately avoid.
#   [A2] update: we KEEP the repo's sign(PnD) 1-D coordinate step (not the
#        paper's raw-gradient P_n(D) step, which can stall with vanishing
#        steps before the stopping criterion is met).
#   [A3] cost: the per-covariate tau_s / gamma_s regressions use a strict,
#        lightweight SL.glmnet library (the 2p secondary LOO parameters),
#        isolated from the heavy shared glmnet+cforest used for Q / g.
#   [A4] Psi_3: CIs use the paper's natural-scale EIC
#        D*_{Psi3} = (D*_{Psi2} - Psi3 * D*_{Psi1}) / Psi1
#        (algebraically equal to the repo's df_log_rr log-ratio IC times the
#        Jacobian Psi3; the repo's log/exp variant is kept commented below).
#
# This file is standalone: it depends only on SuperLearner + glmnet (already
# loaded by R/meta_learners.R) and the small helper .coerce_numeric_df
# defined there. Sourcing it adds compute_te_vim() and the TMLE_* primitives.
############################################################################


# ==========================================================================
# §1. HELPERS (verbatim from vim.R / fit_para.R)
# ==========================================================================

# bound() — clamp to [lower, upper] (vim.R / fit_para.R).
.tevim_bound <- function(x, bounds) {
  lower <- bounds[[1]]
  upper <- if (length(bounds) > 1) bounds[[2]] else 1 - lower
  pmin(pmax(x, lower), upper)
}

# Local numeric-coercion fallback so this file also runs without the harness
# (R/meta_learners.R defines an identical .coerce_numeric_df).
if (!exists(".coerce_numeric_df", mode = "function")) {
  .coerce_numeric_df <- function(X) {
    X <- as.data.frame(X)
    for (col in names(X)) if (!is.numeric(X[[col]])) X[[col]] <- as.numeric(as.factor(X[[col]]))
    X
  }
}


# ==========================================================================
# §2. NUISANCE SUB-FITS  tau_s, gamma_s   (port of fit_para.R::fit_x)
#
# fit_x(outcome = 'tau', para = 'tau_s', ws = s)  regresses tau   on W_{-s}.
# fit_x(outcome = 'tau', para = 'gamma_s', ws = s) regresses tau^2 on W_{-s}.
# Here `s` is the EXCLUDED subset (ws); W_{-s} = setdiff(names(X), s).
# [A3] strict SL.glmnet library; lm/ns and mean fallbacks for the degenerate
#      (<2 column) case, mirroring the project's fit_subcate_default.
# ==========================================================================
.tevim_fit_one <- function(y_train, X_sub, sub_library = "SL.glmnet") {
  X_sub <- as.data.frame(X_sub)
  if (ncol(X_sub) == 0) return(rep(mean(y_train), nrow(X_sub)))

  # glmnet needs >= 2 predictor columns; fall back to a natural-spline lm.
  if (ncol(X_sub) == 1) {
    df_tr <- data.frame(.y = y_train, x = X_sub[[1]])
    deg   <- min(4, length(unique(df_tr$x)) - 1)
    fit   <- tryCatch(
      if (deg >= 3) lm(.y ~ splines::ns(x, df = deg), data = df_tr) else lm(.y ~ x, data = df_tr),
      error = function(e) NULL)
    if (is.null(fit)) return(rep(mean(y_train), nrow(X_sub)))
    pred <- tryCatch(as.numeric(predict(fit, df_tr)), error = function(e) NULL)
    if (is.null(pred)) return(rep(mean(y_train), nrow(X_sub)))
    pred[!is.finite(pred)] <- mean(y_train)
    return(pred)
  }

  fit <- tryCatch(
    SuperLearner::SuperLearner(Y = y_train, X = X_sub, newX = X_sub,
                               family = gaussian(), SL.library = sub_library,
                               verbose = FALSE),
    error = function(e) NULL)
  if (!is.null(fit)) {
    pred <- tryCatch(as.numeric(predict(fit, X_sub, onlySL = TRUE)$pred),
                     error = function(e) NULL)
    if (!is.null(pred) && all(is.finite(pred))) return(pred)
  }
  # Direct cv.glmnet, then mean.
  fit2 <- tryCatch(glmnet::cv.glmnet(as.matrix(X_sub), y_train, alpha = 1, nfolds = 5),
                   error = function(e) NULL)
  if (!is.null(fit2)) {
    pred <- as.numeric(predict(fit2, newx = as.matrix(X_sub), s = "lambda.min"))
    pred[!is.finite(pred)] <- mean(y_train)
    return(pred)
  }
  rep(mean(y_train), nrow(X_sub))
}

# Returns tau_s = E[tau | W_{-s}] and gamma_s = E[tau^2 | W_{-s}].
.tevim_fit_tau_s_gamma_s <- function(X, tau, s, sub_library = "SL.glmnet") {
  X      <- .coerce_numeric_df(X)
  keep   <- setdiff(names(X), s)                 # W_{-s}
  X_sub  <- X[, keep, drop = FALSE]
  tau_s   <- .tevim_fit_one(as.numeric(tau),   X_sub, sub_library)
  gamma_s <- .tevim_fit_one(as.numeric(tau)^2, X_sub, sub_library)
  list(tau_s = tau_s, gamma_s = gamma_s)
}


# ==========================================================================
# §3. TMLE PRIMITIVES (verbatim ports of vim.R, with [A1]/[A2]/[A4])
# ==========================================================================

# --------------------------------------------------------------------------
# TMLE_VIM — Psi_2 = var(tau) - var(tau_s)   (port of vim.R::TMLE_VIM, L675)
# `data` columns: A, Y, pi_hat, mu0_hat, mu1_hat, mua_hat, tau, tau_s, gamma_s.
# --------------------------------------------------------------------------
TMLE_VIM <- function(data, max_it = 1e4, lr = 1e-4, verbose = FALSE) {
  N  <- nrow(data)
  A  <- data$A
  Y  <- data$Y

  QA_0 <- data$mua_hat
  Q1_0 <- data$mu1_hat
  Q0_0 <- data$mu0_hat
  gn   <- data$pi_hat

  tau_0     <- data$tau          ## [A1] = mu1_hat - mu0_hat (S-learner)
  tau_s_0   <- data$tau_s
  gamma_s_0 <- data$gamma_s

  # sig1, sig2 — scale of the two EIC components (for the stopping rule).
  sig1 <- sd(2 * (tau_0 - tau_s_0) * (A / gn - (1 - A) / (1 - gn)) * (Y - QA_0))
  sig2 <- sd(2 * tau_s_0 * (tau_0 - tau_s_0))
  # Guard degenerate (homogeneous) cases so the loop still terminates sanely.
  if (!is.finite(sig1) || sig1 <= 0) sig1 <- Inf
  if (!is.finite(sig2) || sig2 <= 0) sig2 <- Inf

  eps1 <- lr; eps2 <- lr

  i      <- 1
  QA_i   <- QA_0; Q1_i <- Q1_0; Q0_i <- Q0_0
  tau_i  <- tau_0; tau_s_i <- tau_s_0

  while (i <= max_it) {
    # step 1. update Q and tau
    if (i == 1) {
      HQ_1i <- 2 * (tau_i - tau_s_i)
      H1_1i <- 2 * (tau_i - tau_s_i) * (1 / gn)
      H0_1i <- 2 * (tau_i - tau_s_i) * (-1 / (1 - gn))
      HA_1i <- ifelse(A, H1_1i, H0_1i)
      D_1i  <- HA_1i * (Y - QA_i)
      PnD_1i <- mean(D_1i)
    }

    Q1_i <- Q1_i + eps1 * (HQ_1i) * sign(PnD_1i)       ## [A2] sign() coordinate step
    Q0_i <- Q0_i + eps1 * (-HQ_1i) * sign(PnD_1i)
    QA_i <- ifelse(A, Q1_i, Q0_i)
    tau_i <- Q1_i - Q0_i                               ## [A1] tau derived from Q

    # step 2. update tau_s
    H_2i  <- 2 * tau_s_i
    D_2i  <- H_2i * (tau_i - tau_s_i)
    PnD_2i <- mean(D_2i)
    tau_s_i <- tau_s_i + eps2 * H_2i * sign(PnD_2i)

    # recompute D1, D2 and check criteria
    HQ_1i <- 2 * (tau_i - tau_s_i)
    H1_1i <- 2 * (tau_i - tau_s_i) * (1 / gn)
    H0_1i <- 2 * (tau_i - tau_s_i) * (-1 / (1 - gn))
    HA_1i <- ifelse(A, H1_1i, H0_1i)
    D_1i  <- HA_1i * (Y - QA_i)
    PnD_1i <- mean(D_1i)

    H_2i  <- 2 * tau_s_i
    D_2i  <- H_2i * (tau_i - tau_s_i)
    PnD_2i <- mean(D_2i)

    c1 <- abs(PnD_1i) <= sig1 / (sqrt(N) * log(N))
    c2 <- abs(PnD_2i) <= sig2 / (sqrt(N) * log(N))
    if (c1 && c2) break
    i <- i + 1
  }

  tau_star   <- tau_i
  tau_s_star <- tau_s_i
  if (i >= max_it && verbose) warning("Max iterations reached in TMLE_VIM")
  if (verbose) message(sprintf("TMLE_VIM steps: %d", i))

  # update gamma_s by a one-parameter linear (intercept) fluctuation:
  # gamma_s* = gamma_s_0 + mean(tau_star^2 - gamma_s_0)  =>  mean(gamma_s*) = mean(tau_star^2)
  suppressWarnings({
    linearUpdate <- glm(tau_star^2 ~ 1, offset = gamma_s_0, family = "gaussian")
  })
  eps3 <- coef(linearUpdate)
  gamma_s_star <- gamma_s_0 + eps3

  theta_s_star <- mean(gamma_s_star - tau_s_star^2)     # Psi_2 (targeted)

  # Influence curve at the INITIAL fits (repo convention "use initial est for ic").
  theta_s_0 <- mean(gamma_s_0 - tau_s_0^2)
  ic <- 2 * (tau_0 - tau_s_0) * (A / gn - (1 - A) / (1 - gn)) * (Y - QA_0) +
        (tau_0 - tau_s_0)^2 - theta_s_0
  se <- sqrt(var(ic) / N)

  list(coef = theta_s_star, std_err = se,
       ci_l = theta_s_star - 1.96 * se, ci_u = theta_s_star + 1.96 * se,
       ic = ic, n_iter = i)
}

# --------------------------------------------------------------------------
# TMLE_VTE — Psi_1 = var(tau)   (port of vim.R::TMLE_VTE, L1009)
# `data` columns: A, Y, pi_hat, mu0_hat, mu1_hat, mua_hat, tau.
# --------------------------------------------------------------------------
TMLE_VTE <- function(data, max_it = 1e4, lr = 1e-4, verbose = FALSE) {
  N  <- nrow(data)
  A  <- data$A
  Y  <- data$Y

  QA_0 <- data$mua_hat
  Q1_0 <- data$mu1_hat
  Q0_0 <- data$mu0_hat
  gn   <- data$pi_hat

  tau_0  <- data$tau                                    ## [A1] S-learner
  ate_0  <- rep(mean(tau_0), N)

  sig1 <- sd(2 * (tau_0 - ate_0) * (A / gn - (1 - A) / (1 - gn)) * (Y - QA_0))
  if (!is.finite(sig1) || sig1 <= 0) sig1 <- Inf
  eps1 <- lr

  i     <- 1
  QA_i  <- QA_0; Q1_i <- Q1_0; Q0_i <- Q0_0
  tau_i <- tau_0; ate_i <- ate_0

  while (i <= max_it) {
    if (i == 1) {
      HQ_1i <- 2 * (tau_i - ate_i)
      H1_1i <- 2 * (tau_i - ate_i) * (1 / gn)
      H0_1i <- 2 * (tau_i - ate_i) * (-1 / (1 - gn))
      HA_1i <- ifelse(A, H1_1i, H0_1i)
      D_1i  <- HA_1i * (Y - QA_i)
      PnD_1i <- mean(D_1i)
    }

    Q1_i <- Q1_i + eps1 * (HQ_1i) * sign(PnD_1i)        ## [A2]
    Q0_i <- Q0_i + eps1 * (-HQ_1i) * sign(PnD_1i)
    QA_i <- ifelse(A, Q1_i, Q0_i)
    tau_i <- Q1_i - Q0_i
    ate_i <- rep(mean(tau_i), N)

    HQ_1i <- 2 * (tau_i - ate_i)
    H1_1i <- 2 * (tau_i - ate_i) * (1 / gn)
    H0_1i <- 2 * (tau_i - ate_i) * (-1 / (1 - gn))
    HA_1i <- ifelse(A, H1_1i, H0_1i)
    D_1i  <- HA_1i * (Y - QA_i)
    PnD_1i <- mean(D_1i)

    c1 <- abs(PnD_1i) <= sig1 / (sqrt(N) * log(N))
    if (c1) break
    i <- i + 1
  }

  tau_star <- tau_i
  ate_star <- rep(mean(tau_i), N)
  if (i >= max_it && verbose) warning("Max iterations reached in TMLE_VTE")

  theta_s_star <- mean((tau_star - ate_star)^2)         # Psi_1 (targeted VTE)
  theta_s_0    <- mean((tau_0 - ate_0)^2)
  ic <- 2 * (tau_0 - ate_0) * (A / gn - (1 - A) / (1 - gn)) * (Y - QA_0) +
        (tau_0 - ate_0)^2 - theta_s_0
  se <- sqrt(var(ic) / N)

  list(coef = theta_s_star, std_err = se,
       ci_l = theta_s_star - 1.96 * se, ci_u = theta_s_star + 1.96 * se,
       ic = ic, n_iter = i)
}

# --------------------------------------------------------------------------
# compute_psi3 — Psi_3 = Psi_2 / Psi_1   (scaled VIMb).
# [A4] paper natural-scale EIC:  D*_{Psi3} = (ic_{Psi2} - Psi3 * ic_{Psi1}) / Psi1.
#      (Repo TMLE_VIM2 instead builds the log-ratio IC ic_{Psi2}/Psi2 -
#       ic_{Psi1}/Psi1 and exponentiates; equal up to the Jacobian Psi3.)
# --------------------------------------------------------------------------
compute_psi3 <- function(res_vim, res_vte) {
  psi2 <- res_vim$coef   # Psi_2 (VIMa)
  psi1 <- res_vte$coef   # Psi_1 (VTE)
  N    <- length(res_vim$ic)
  if (!is.finite(psi1) || abs(psi1) < 1e-12) {
    return(list(coef = NA_real_, std_err = NA_real_, ci_l = NA_real_, ci_u = NA_real_))
  }
  psi3 <- psi2 / psi1
  ic3  <- (res_vim$ic - psi3 * res_vte$ic) / psi1
  se   <- sqrt(var(ic3) / N)

  # # Repo log-ratio variant (kept for reference):
  # ic_log <- res_vim$ic / psi2 - res_vte$ic / psi1
  # se_log <- sqrt(var(ic_log) / N)
  # ci on natural scale = exp(log(psi3) +/- 1.96 * se_log)

  list(coef = psi3, std_err = se,
       ci_l = psi3 - 1.96 * se, ci_u = psi3 + 1.96 * se, ic = ic3)
}


# ==========================================================================
# §4. DRIVER:  compute_te_vim()
#
# Builds the S-learner CATE from the SHARED nuisances [A1], computes the VTE
# (Psi_1) once, then loops the requested subsets (default = each single
# covariate, "leave-one-out") computing Psi_2 (VIMa) and Psi_3 (VIMb) per
# covariate. Returns one tidy row per subset for the harness `tevims_scores`.
# ==========================================================================
compute_te_vim <- function(X, Y, trt, nuis,
                           subsets = NULL,          # list of excluded sets; NULL = LOO each covar
                           sub_library = "SL.glmnet",  # [A3]
                           max_it = 1e4, lr = 1e-4,
                           verbose = FALSE) {
  X   <- .coerce_numeric_df(X)
  Y   <- as.numeric(Y); trt <- as.numeric(trt); N <- nrow(X)
  covs <- names(X)
  if (is.null(subsets)) subsets <- as.list(covs)   # LOO: one covariate per subset

  # [A1] S-learner CATE from the shared nuisances (tau = Q1 - Q0).
  tau0 <- as.numeric(nuis$mu1_hat - nuis$mu0_hat)

  # te_vim bounds g to [0.025, 0.975]; rebind for fidelity (harness clips [0.01,0.99]).
  gn <- .tevim_bound(as.numeric(nuis$e_hat), c(0.025, 0.975))
  QA <- ifelse(trt == 1, nuis$mu1_hat, nuis$mu0_hat)

  base <- data.frame(A = trt, Y = Y, pi_hat = gn,
                     mu0_hat = as.numeric(nuis$mu0_hat),
                     mu1_hat = as.numeric(nuis$mu1_hat),
                     mua_hat = as.numeric(QA), tau = tau0)

  # VTE (Psi_1) is independent of s — compute ONCE.
  res_vte <- TMLE_VTE(base, max_it = max_it, lr = lr, verbose = verbose)

  rows <- list()
  for (s in subsets) {
    fit  <- .tevim_fit_tau_s_gamma_s(X, tau0, s, sub_library = sub_library)
    df_s <- base
    df_s$tau_s   <- fit$tau_s
    df_s$gamma_s <- fit$gamma_s

    res_vim  <- TMLE_VIM(df_s, max_it = max_it, lr = lr, verbose = verbose)  # Psi_2
    res_psi3 <- compute_psi3(res_vim, res_vte)                               # Psi_3

    rows[[length(rows) + 1]] <- data.frame(
      variable  = paste(s, collapse = "+"),
      mode      = "LOO",
      estimator = "TMLE",
      vima      = res_vim$coef,  vima_se   = res_vim$std_err,
      vima_ci_l = res_vim$ci_l,  vima_ci_u = res_vim$ci_u,
      vimb      = res_psi3$coef, vimb_se   = res_psi3$std_err,
      vimb_ci_l = res_psi3$ci_l, vimb_ci_u = res_psi3$ci_u,
      vte       = res_vte$coef,  vte_se    = res_vte$std_err,
      n_iter_vim = res_vim$n_iter, n_iter_vte = res_vte$n_iter,
      stringsAsFactors = FALSE)
  }
  if (length(rows) == 0) return(NULL)
  out <- do.call(rbind, rows); rownames(out) <- NULL
  out
}
