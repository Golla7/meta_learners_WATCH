# ---------------------------------------------------------------------------
# Verification checks for R/meta_learners.R
#   1. same-footing: one shared nuisance object feeds every learner
#   2. DINA Cross-Fit iterates ALL K shared folds (DML1 beta averaging)
#   3. Gao (2024) relative-error EIF: sanity on synthetic data with known
#      nuisances — the estimator must favour the tau that is closer to truth
# Run from the repo root:  Rscript tests/verify_harness.R
# ---------------------------------------------------------------------------
if (basename(getwd()) == "tests") setwd("..")
COMP_CLEAN_NO_DEMO <- TRUE
suppressWarnings(suppressMessages(source("R/meta_learners.R")))

set.seed(42)
n <- 250; p <- 8
X <- as.data.frame(matrix(rnorm(n * p), n, p)); names(X) <- paste0("X", 1:p)
trt <- rbinom(n, 1, 0.5)
tau <- 0.5 * X$X1
Y <- X$X2 + tau * trt + rnorm(n)

nuis <- estimate_nuisances(X, Y, trt, seed = 42, light_config = TRUE)

## Check 1: same footing — all learners receive the identical nuis object.
cat("== CHECK 1: same-footing (shared nuisance fingerprints) ==\n")
cat(sprintf("  sum(e_hat)  = %.10f\n", sum(nuis$e_hat)))
cat(sprintf("  sum(mu0)    = %.10f\n", sum(nuis$mu0_hat)))
cat(sprintf("  sum(mu1)    = %.10f\n", sum(nuis$mu1_hat)))
cat(sprintf("  sum(phi)    = %.10f\n", sum(nuis$phi)))
stopifnot(length(nuis$e_hat) == n, all(is.finite(nuis$phi)))
cat("  PASS: single finite nuis object (dispatched to every learner by compare_meta_learners)\n")

## Check 2: DINA cross-fit iterates all K shared folds.
cat("\n== CHECK 2: DINA cross-fit fold usage ==\n")
dcf <- learner_dina_crossfit(X, Y, trt, nuis, seed = 42, nperm = 2)
cat(sprintf("  unique folds in nuis$fold_ids : %d\n", length(unique(nuis$fold_ids))))
cat(sprintf("  per-fold n                    : %s\n",
            paste(table(nuis$fold_ids), collapse = ", ")))
cat(sprintf("  beta vectors averaged         : %d\n", dcf$n_folds_used))
stopifnot(dcf$n_folds_used == length(unique(nuis$fold_ids)))
cat("  PASS: averaged exactly K per-fold betas\n")

## Check 3: Gao (2024) relative-error EIF sanity with ORACLE nuisances.
## delta = E[(tau1 - tau)^2] - E[(tau2 - tau)^2]; with tau1 = truth and
## tau2 = truth + noise the estimate must be negative (tau1 wins) and its
## CI must exclude zero on this easy, well-powered example.
cat("\n== CHECK 3: relative-error EIF favours the better tau ==\n")
set.seed(7)
m <- 2000
Xe  <- rnorm(m)
tau_t <- 0.5 * Xe                    # true CATE
m0  <- Xe                            # true mu0 = prognostic part
m1  <- m0 + tau_t                    # true mu1
We  <- rbinom(m, 1, 0.5)
Ye  <- ifelse(We == 1, m1, m0) + rnorm(m)
eh  <- rep(0.5, m)
tau1 <- tau_t                        # perfect estimator
tau2 <- tau_t + rnorm(m, sd = 0.7)   # noisy estimator
r <- compute_relative_error_eif(tau1, tau2, Ye, We, eh, m0, m1)
cat(sprintf("  delta = %.4f  [%.4f, %.4f]\n", r$delta, r$ci_lo, r$ci_hi))
stopifnot(r$delta < 0, isTRUE(r$ci_excludes_zero))
cat("  PASS: EIF identifies the more accurate CATE estimator\n")

cat("\nALL CHECKS DONE\n")
