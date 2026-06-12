#!/usr/bin/env Rscript
# tests/test_te_vim.R — standalone verification of the te_vim TMLE port.
# Run from the repo root:  Rscript tests/test_te_vim.R
#
# (1) Reproduces the paper's section 4.1 DGP and runs compute_te_vim() with
#     ORACLE nuisances (so the TMLE targeting is tested in isolation).
# (2) Cross-checks TMLE_VIM() against an INLINE verbatim copy of the original
#     HaodongL/te_vim TMLE_VIM on identical df_fit — they must agree to ~1e-10.

if (basename(getwd()) == "tests") setwd("..")
suppressPackageStartupMessages({ library(SuperLearner); library(glmnet) })
source("R/te_vim.R")

set.seed(20240609)

# ---- Paper section 4.1 DGP ------------------------------------------------
expit <- function(z) 1 / (1 + exp(-z))
gen <- function(n) {
  W1 <- runif(n, -1, 1); W2 <- runif(n, -1, 1)
  p  <- expit(0.1 * W1 * W2 - 0.4 * W1)
  A  <- rbinom(n, 1, p)
  tau <- W1^2 * (W1 + 7/5) + (5 * W2 / 3)^2
  prog <- W1 * W2 + 2 * W2^2 - W1
  muY <- A * tau + prog
  Y  <- rnorm(n, muY, 1)
  list(X = data.frame(W1 = W1, W2 = W2), Y = Y, trt = A,
       tau = tau, e = p, mu1 = tau + prog, mu0 = prog)
}

d <- gen(2000)
# Oracle nuisances (S-learner tau0 = mu1 - mu0 = true tau).
nuis <- list(e_hat = d$e, mu0_hat = d$mu0, mu1_hat = d$mu1)

cat("===== (1) compute_te_vim with ORACLE nuisances (n=2000) =====\n")
res <- compute_te_vim(d$X, d$Y, d$trt, nuis,
                      subsets = list("W1", "W2"),
                      sub_library = "SL.glmnet", verbose = TRUE)
print(res[, c("variable","vima","vima_ci_l","vima_ci_u",
              "vimb","vte","n_iter_vim","n_iter_vte")])

stopifnot(all(is.finite(res$vima)),
          all(res$n_iter_vim < 1e4),     # stopping criterion actually trips
          res$n_iter_vte[1] < 1e4,
          all(res$vte > 0))
cat("\n[OK] VIMa finite, VTE > 0, TMLE loops converged (i < max_it).\n")
cat(sprintf("  Ranking by VIMa:  %s\n",
            paste(res$variable[order(-res$vima)], collapse = " > ")))

# ---- (2) Faithfulness cross-check vs the ORIGINAL repo TMLE_VIM ------------
# Verbatim copy of HaodongL/te_vim R/est_function/vim.R::TMLE_VIM (master).
TMLE_VIM_REPO <- function(data, max_it = 600, lr = 1e-4, logtrans = FALSE){
  N <- nrow(data); A <- data$A; Y <- data$Y; po <- data$po
  QA_0 <- data$mua_hat; Q1_0 <- data$mu1_hat; Q0_0 <- data$mu0_hat; gn <- data$pi_hat
  tau_0 <- data$tau; tau_s_0 <- data$tau_s; gamma_s_0 <- data$gamma_s
  sig1 <- sd(2*(tau_0 - tau_s_0)*(A/gn - (1-A)/(1-gn))*(Y - QA_0))
  sig2 <- sd(2*tau_s_0*(tau_0 - tau_s_0))
  eps1 <- lr; eps2 <- lr
  i <- 1; QA_i <- QA_0; Q1_i <- Q1_0; Q0_i <- Q0_0; tau_i <- tau_0; tau_s_i <- tau_s_0
  while(i <= max_it){
    if (i == 1){
      HQ_1i <- 2*(tau_i - tau_s_i)
      H1_1i <- 2*(tau_i - tau_s_i)*(1/gn); H0_1i <- 2*(tau_i - tau_s_i)*(-1/(1-gn))
      HA_1i <- ifelse(A, H1_1i, H0_1i); D_1i <- HA_1i*(Y - QA_i); PnD_1i <- mean(D_1i)
    }
    Q1_i <- Q1_i + eps1*(HQ_1i)*sign(PnD_1i); Q0_i <- Q0_i + eps1*(-HQ_1i)*sign(PnD_1i)
    QA_i <- ifelse(A, Q1_i, Q0_i); tau_i <- Q1_i - Q0_i
    H_2i <- 2*tau_s_i; D_2i <- H_2i*(tau_i - tau_s_i); PnD_2i <- mean(D_2i)
    tau_s_i <- tau_s_i + eps2*H_2i*sign(PnD_2i)
    HQ_1i <- 2*(tau_i - tau_s_i)
    H1_1i <- 2*(tau_i - tau_s_i)*(1/gn); H0_1i <- 2*(tau_i - tau_s_i)*(-1/(1-gn))
    HA_1i <- ifelse(A, H1_1i, H0_1i); D_1i <- HA_1i*(Y - QA_i); PnD_1i <- mean(D_1i)
    H_2i <- 2*tau_s_i; D_2i <- H_2i*(tau_i - tau_s_i); PnD_2i <- mean(D_2i)
    c1 <- abs(PnD_1i) <= sig1/(sqrt(N)*log(N)); c2 <- abs(PnD_2i) <= sig2/(sqrt(N)*log(N))
    if (c1 & c2){ break }
    i <- i + 1
  }
  QA_star <- QA_i; tau_star <- tau_i; tau_s_star <- tau_s_i
  suppressWarnings({ linearUpdate <- glm(tau_star^2 ~ 1, offset = gamma_s_0, family = "gaussian") })
  eps3 <- coef(linearUpdate); gamma_s_star <- gamma_s_0 + eps3
  theta_s_star <- mean(gamma_s_star - tau_s_star^2)
  theta_s_0 <- mean(gamma_s_0 - tau_s_0^2)
  ic <- 2*(tau_0 - tau_s_0)*(A/gn - (1-A)/(1-gn))*(Y-QA_0) + (tau_0 - tau_s_0)^2 - theta_s_0
  se <- sqrt(var(ic)/N)
  list(coef = theta_s_star, std_err = se, ci_l = theta_s_star - 1.96*se,
       ci_u = theta_s_star + 1.96*se, ic = ic)
}

cat("\n===== (2) Port vs original repo TMLE_VIM on identical df_fit =====\n")
# Build one df_fit exactly as the driver does for s = {W1}.
tau0 <- nuis$mu1_hat - nuis$mu0_hat
gn   <- pmax(0.025, pmin(nuis$e_hat, 0.975))
QA   <- ifelse(d$trt == 1, nuis$mu1_hat, nuis$mu0_hat)
ff   <- .tevim_fit_tau_s_gamma_s(d$X, tau0, "W1", "SL.glmnet")
df_fit <- data.frame(A = d$trt, Y = d$Y, pi_hat = gn,
                     mu0_hat = nuis$mu0_hat, mu1_hat = nuis$mu1_hat, mua_hat = QA,
                     po = 0, tau = tau0, tau_s = ff$tau_s, gamma_s = ff$gamma_s)

mine <- TMLE_VIM(df_fit, max_it = 1e4, lr = 1e-4)
repo <- TMLE_VIM_REPO(df_fit, max_it = 1e4, lr = 1e-4)
cat(sprintf("  port  Psi_2 = %.10f  (se %.10f)\n", mine$coef, mine$std_err))
cat(sprintf("  repo  Psi_2 = %.10f  (se %.10f)\n", repo$coef, repo$std_err))
stopifnot(abs(mine$coef - repo$coef)       < 1e-9,
          abs(mine$std_err - repo$std_err) < 1e-9)
cat("[OK] Port reproduces the original TMLE_VIM to < 1e-9.\n")

cat("\nALL TE-VIM PORT TESTS PASSED.\n")
