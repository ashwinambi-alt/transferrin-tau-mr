# =============================================================================
# TRANSFERRIN AND TAU BIOMARKERS: MULTIVARIABLE MENDELIAN RANDOMIZATION
# =============================================================================
# Reproducibility test: 15-run stability validation of primary results
# Paper: "Iron transport capacity causally reduces circulating and CSF tau
#         biomarkers independent of serum iron: a multivariable Mendelian
#         randomization study"
# Author: Ashwin Ambi
# Preprint: [TO BE ADDED ON SUBMISSION]
# Journal submission: JNNP (target)
#
# R version: 4.5.3 (2026-03-11 ucrt)
# Key packages: TwoSampleMR 0.7.1, MendelianRandomization, dplyr, ggplot2
#
# All GWAS summary data accessed via IEU OpenGWAS API (https://opengwas.io/)
# and GWAS Catalog (https://www.ebi.ac.uk/gwas/)
# No individual-level data used.
# =============================================================================

##############################################################################
# REPRODUCIBILITY_TEST.R — Verify results are stable across 15 runs
#
# MR methods (IVW, Egger, WM, WMode) are DETERMINISTIC — they produce
# identical results every time given the same input data.
# MR-PRESSO uses permutations — varies by seed.
#
# This script:
#   1. Loads the harmonised data from the master pipeline
#   2. Runs the headline analysis (transferrin → tau) 15 times
#   3. Records IVW beta, SE, p-value each run
#   4. Records MR-PRESSO global p each run (with different seeds)
#   5. Runs MVMR 15 times
#   6. Prints summary confirming reproducibility
##############################################################################

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(MendelianRandomization)
  library(MRPRESSO)
})

DATA_DIR <- normalizePath(file.path(getwd(), "..", "data"), mustWork = TRUE)
RESULTS_DIR <- normalizePath(file.path(getwd(), "..", "results"), mustWork = TRUE)

cat("==============================================================\n")
cat("  REPRODUCIBILITY TEST — 15 runs\n")
cat("  Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("==============================================================\n\n")

# Load data from master pipeline
load(file.path(DATA_DIR, "FINAL_analysis_results.RData"))

dat_tt <- all_harmonised[["transferrin_circ_total_tau"]]
if (is.null(dat_tt)) stop("Harmonised data for transferrin -> tau not found. Run MASTER_PIPELINE.R first.")

N_RUNS <- 15

MR_METHODS <- c("mr_ivw", "mr_egger_regression", "mr_weighted_median",
                "mr_weighted_mode", "mr_simple_mode")

# ============================================================
# 1. Transferrin → Tau: IVW + all methods (15 runs)
# ============================================================

cat("=== Test 1: Transferrin -> Tau MR (deterministic methods) ===\n\n")

ivw_results <- data.frame(run = integer(), beta = numeric(), se = numeric(),
                           pval = numeric(), nsnp = integer())
all_method_results <- list()

for (i in 1:N_RUNS) {
  set.seed(i * 7)  # Different seed each run (shouldn't matter for deterministic methods)
  mr_res <- mr(dat_tt, method_list = MR_METHODS)
  ivw <- mr_res[mr_res$method == "Inverse variance weighted", ]
  ivw_results <- rbind(ivw_results, data.frame(
    run = i, beta = ivw$b, se = ivw$se, pval = ivw$pval, nsnp = ivw$nsnp))
  all_method_results[[i]] <- mr_res
}

cat("IVW results across", N_RUNS, "runs:\n")
cat(sprintf("  %-5s %12s %12s %15s %6s\n", "Run", "Beta", "SE", "P-value", "SNPs"))
cat(paste(rep("-", 55), collapse = ""), "\n")
for (i in 1:N_RUNS) {
  cat(sprintf("  %-5d %12.6f %12.6f %15s %6d\n",
    ivw_results$run[i], ivw_results$beta[i], ivw_results$se[i],
    formatC(ivw_results$pval[i], format = "e", digits = 4), ivw_results$nsnp[i]))
}

# Check consistency
beta_range <- diff(range(ivw_results$beta))
se_range <- diff(range(ivw_results$se))
p_range <- diff(range(ivw_results$pval))

cat("\nConsistency check:\n")
cat("  Beta range:", beta_range, ifelse(beta_range == 0, " IDENTICAL", " DIFFERS"), "\n")
cat("  SE range:  ", se_range, ifelse(se_range == 0, " IDENTICAL", " DIFFERS"), "\n")
cat("  P range:   ", p_range, ifelse(p_range == 0, " IDENTICAL", " DIFFERS"), "\n")

if (beta_range == 0 & se_range == 0 & p_range == 0) {
  cat("\n  PASS: IVW results are perfectly reproducible (deterministic).\n")
} else {
  cat("\n  FLAG: Results vary — investigate.\n")
}

# Check all 5 methods
cat("\nAll 5 methods (Run 1 vs Run 15):\n")
for (m in MR_METHODS) {
  r1 <- all_method_results[[1]][all_method_results[[1]]$method == mr_res$method[mr_res$method == gsub("mr_", "", m) | TRUE][1], ]
  # Simpler: just compare
}
m1 <- all_method_results[[1]]
m15 <- all_method_results[[15]]
for (j in 1:nrow(m1)) {
  match <- m1$b[j] == m15$b[j] & m1$se[j] == m15$se[j]
  cat(sprintf("  %-30s Run1: %.6f  Run15: %.6f  %s\n",
    m1$method[j], m1$b[j], m15$b[j], ifelse(match, "IDENTICAL", "DIFFERS")))
}

# ============================================================
# 2. MR-PRESSO (stochastic — varies by seed)
# ============================================================

cat("\n\n=== Test 2: MR-PRESSO (stochastic, 15 seeds) ===\n\n")

presso_results <- data.frame(run = integer(), seed = integer(),
                              global_p = numeric(), estimate = numeric())

for (i in 1:N_RUNS) {
  seed_val <- i * 13 + 7
  set.seed(seed_val)
  tryCatch({
    res <- mr_presso(
      BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
      SdOutcome = "se.outcome", SdExposure = "se.exposure",
      OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
      data = dat_tt, NbDistribution = 10000, SignifThreshold = 0.05)
    gp <- res[["MR-PRESSO results"]][["Global Test"]]$Pvalue
    est <- res[["Main MR results"]][["Causal Estimate"]][1]
    presso_results <- rbind(presso_results, data.frame(
      run = i, seed = seed_val, global_p = gp, estimate = est))
    cat(sprintf("  Run %2d (seed=%4d): global p = %.4f, estimate = %.5f\n", i, seed_val, gp, est))
  }, error = function(e) {
    cat(sprintf("  Run %2d: ERROR\n", i))
  })
}

if (nrow(presso_results) > 0) {
  cat("\nMR-PRESSO summary:\n")
  cat("  Global p: mean =", round(mean(presso_results$global_p, na.rm = TRUE), 4),
      ", min =", round(min(presso_results$global_p, na.rm = TRUE), 4),
      ", max =", round(max(presso_results$global_p, na.rm = TRUE), 4), "\n")
  cat("  Estimate: mean =", round(mean(presso_results$estimate, na.rm = TRUE), 5),
      ", sd =", round(sd(presso_results$estimate, na.rm = TRUE), 6), "\n")
  cat("  All runs global p > 0.05?", all(presso_results$global_p > 0.05, na.rm = TRUE), "\n")
  cat("  Verdict:", ifelse(all(presso_results$global_p > 0.05, na.rm = TRUE),
    "PASS — no outliers detected in any run", "CHECK — some runs flagged outliers"), "\n")
}

# ============================================================
# 3. MVMR (deterministic — 15 runs)
# ============================================================

cat("\n\n=== Test 3: MVMR (15 runs) ===\n\n")

# Reconstruct MVMR input from saved data
tf_inst <- all_instruments[["transferrin"]]
fe_inst <- all_instruments[["iron"]]

mvmr_tf_results <- data.frame(run = integer(), beta = numeric(),
                               se = numeric(), pval = numeric())

if (!is.null(tf_inst) && !is.null(fe_inst)) {
  # We need the raw MVMR matrices — reconstruct from the FINAL data
  # Load MVMR result to confirm we have the input
  Sys.setenv(OPENGWAS_JWT = "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJhYW1iaUBzYWludHBldGVycy5lZHUiLCJpYXQiOjE3NzQ2NzY0MzcsImV4cCI6MTc3NTg4NjAzN30.bAM7zgV0nx3DU8FlLtK7z5n9lpPpZUZOqeTwoK2417Al0vbGluvG57Io55Q6Xe-3zcbmD9yea4s8f-aK9A7P7RuLXCZM_wn3KoZFP1FFf5cx12tPFkVc1BVCoDJu2TxSPrV3npMEbzr15IjXc1Wb9tTvDrdELgYFjk0261QU_YUAmhhxFMqb3NqazeDTQleiK5eLpclw2z_sGo4_m8adB1apARQHsIClYsZuY263b38PIEIaY6bfBNuh0Txz5hkrBppUo2IkQCf7zGgHLF4c98ZJmsNyNEWVnGpZ9z_JZT_eff-fd5NV8lRy6TbVkZL2v06E43kjAUf0MWCTri5QCQ")
  library(ieugwasr)

  all_snps <- unique(c(tf_inst$SNP, fe_inst$SNP))
  tf_exp <- tryCatch(extract_outcome_data(snps = all_snps, outcomes = "ieu-a-1052"), error = function(e) NULL)
  fe_exp <- tryCatch(extract_outcome_data(snps = all_snps, outcomes = "ieu-a-1049"), error = function(e) NULL)
  tau_out <- tryCatch(extract_outcome_data(snps = all_snps, outcomes = "ebi-a-GCST90095138"), error = function(e) NULL)

  if (!is.null(tf_exp) && !is.null(fe_exp) && !is.null(tau_out)) {
    common <- Reduce(intersect, list(tf_exp$SNP, fe_exp$SNP, tau_out$SNP))
    bx <- cbind(tf_exp$beta.outcome[match(common, tf_exp$SNP)],
                 fe_exp$beta.outcome[match(common, fe_exp$SNP)])
    bxse <- cbind(tf_exp$se.outcome[match(common, tf_exp$SNP)],
                   fe_exp$se.outcome[match(common, fe_exp$SNP)])
    by_vec <- tau_out$beta.outcome[match(common, tau_out$SNP)]
    byse_vec <- tau_out$se.outcome[match(common, tau_out$SNP)]
    ok <- complete.cases(bx, by_vec)

    if (sum(ok) >= 5) {
      mvmr_input <- mr_mvinput(bx = bx[ok, ], bxse = bxse[ok, ],
        by = by_vec[ok], byse = byse_vec[ok],
        exposure = c("Transferrin", "Serum Iron"))

      for (i in 1:N_RUNS) {
        set.seed(i * 11)
        mvmr_res <- mr_mvivw(mvmr_input)
        mvmr_tf_results <- rbind(mvmr_tf_results, data.frame(
          run = i, beta = mvmr_res@Estimate[1],
          se = mvmr_res@StdError[1], pval = mvmr_res@Pvalue[1]))
      }

      cat("MVMR Transferrin results across", N_RUNS, "runs:\n")
      cat(sprintf("  %-5s %12s %12s %15s\n", "Run", "Beta", "SE", "P-value"))
      cat(paste(rep("-", 48), collapse = ""), "\n")
      for (i in 1:N_RUNS) {
        cat(sprintf("  %-5d %12.6f %12.6f %15s\n",
          mvmr_tf_results$run[i], mvmr_tf_results$beta[i], mvmr_tf_results$se[i],
          formatC(mvmr_tf_results$pval[i], format = "e", digits = 4)))
      }

      mvmr_range <- diff(range(mvmr_tf_results$beta))
      cat("\n  Beta range:", mvmr_range, ifelse(mvmr_range == 0, " IDENTICAL", " DIFFERS"), "\n")
      cat("  PASS:", ifelse(mvmr_range == 0, "MVMR is perfectly reproducible.", "Check variation."), "\n")
    }
  }
}

# ============================================================
# FINAL SUMMARY
# ============================================================

cat("\n\n==============================================================\n")
cat("  REPRODUCIBILITY SUMMARY\n")
cat("==============================================================\n\n")

cat("Test 1 — IVW (deterministic):\n")
cat("  Beta:", unique(ivw_results$beta), "\n")
cat("  SE:  ", unique(ivw_results$se), "\n")
cat("  P:   ", formatC(unique(ivw_results$pval), format = "e", digits = 4), "\n")
cat("  Verdict:", ifelse(beta_range == 0, "PERFECTLY REPRODUCIBLE", "VARIES"), "\n")

cat("\nTest 2 — MR-PRESSO (stochastic):\n")
if (nrow(presso_results) > 0) {
  cat("  Global p range: [", round(min(presso_results$global_p, na.rm = TRUE), 4),
      ",", round(max(presso_results$global_p, na.rm = TRUE), 4), "]\n")
  cat("  All runs consistent (no outliers)?", all(presso_results$global_p > 0.05, na.rm = TRUE), "\n")
}

cat("\nTest 3 — MVMR (deterministic):\n")
if (nrow(mvmr_tf_results) > 0) {
  cat("  Beta:", unique(mvmr_tf_results$beta), "\n")
  cat("  P:   ", formatC(unique(mvmr_tf_results$pval), format = "e", digits = 4), "\n")
  cat("  Verdict:", ifelse(diff(range(mvmr_tf_results$beta)) == 0, "PERFECTLY REPRODUCIBLE", "VARIES"), "\n")
}

# Save reproducibility report
repro_report <- list(
  ivw = ivw_results,
  presso = presso_results,
  mvmr = mvmr_tf_results,
  timestamp = Sys.time(),
  R_version = R.version.string
)
save(repro_report, file = file.path(DATA_DIR, "reproducibility_report.RData"))
write.csv(ivw_results, file.path(RESULTS_DIR, "primary", "reproducibility_IVW_15runs.csv"), row.names = FALSE)
write.csv(presso_results, file.path(RESULTS_DIR, "primary", "reproducibility_PRESSO_15runs.csv"), row.names = FALSE)
write.csv(mvmr_tf_results, file.path(RESULTS_DIR, "primary", "reproducibility_MVMR_15runs.csv"), row.names = FALSE)

cat("\n\nReproducibility data saved.\n")
cat("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
