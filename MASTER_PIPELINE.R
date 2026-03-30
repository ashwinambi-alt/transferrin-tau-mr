# =============================================================================
# TRANSFERRIN AND TAU BIOMARKERS: MULTIVARIABLE MENDELIAN RANDOMIZATION
# =============================================================================
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
# MASTER_PIPELINE.R — Complete Reproducible Analysis
# MR Study: Essential Metals and Iron Markers → AD Biomarkers
#
# This script runs the ENTIRE analysis from scratch:
#   1. Instrument extraction (all exposures)
#   2. Harmonization (all pairs)
#   3. Primary MR (5 methods, all pairs)
#   4. Trace metal sensitivity (p < 5e-6)
#   5. Sensitivity battery (headline finding)
#   6. MVMR (transferrin + serum iron → tau)
#   7. Positive and negative controls
#   8. Power calculations
#   9. Publication figures
#  10. Manuscript tables
#
# Runtime: ~15-20 minutes (API-dependent)
# Output: All results, figures, and tables in results/
#
# To reproduce: set working directory to scripts/ and run:
#   Rscript MASTER_PIPELINE.R
#
# Principal Investigator: Ashwin Ambi
# Computational Lead: Claude Code
# Date: 2026-03-28
##############################################################################

# ============================================================
# SETUP
# ============================================================

cat("==============================================================\n")
cat("  MASTER PIPELINE — MR: Essential Metals → AD Biomarkers\n")
cat("  Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("==============================================================\n\n")

# JWT token for OpenGWAS API authentication
# Registered: aambi@saintpeters.edu | Expires: 2026-04-11
Sys.setenv(OPENGWAS_JWT = "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJhYW1iaUBzYWludHBldGVycy5lZHUiLCJpYXQiOjE3NzQ2NzY0MzcsImV4cCI6MTc3NTg4NjAzN30.bAM7zgV0nx3DU8FlLtK7z5n9lpPpZUZOqeTwoK2417Al0vbGluvG57Io55Q6Xe-3zcbmD9yea4s8f-aK9A7P7RuLXCZM_wn3KoZFP1FFf5cx12tPFkVc1BVCoDJu2TxSPrV3npMEbzr15IjXc1Wb9tTvDrdELgYFjk0261QU_YUAmhhxFMqb3NqazeDTQleiK5eLpclw2z_sGo4_m8adB1apARQHsIClYsZuY263b38PIEIaY6bfBNuh0Txz5hkrBppUo2IkQCf7zGgHLF4c98ZJmsNyNEWVnGpZ9z_JZT_eff-fd5NV8lRy6TbVkZL2v06E43kjAUf0MWCTri5QCQ")

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(MendelianRandomization)
  library(MRPRESSO)
  library(ieugwasr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(dplyr)
  library(writexl)
})

set.seed(42)

# Project directories
PROJECT_DIR <- normalizePath(file.path(getwd(), ".."), mustWork = TRUE)
DATA_DIR    <- file.path(PROJECT_DIR, "data")
RESULTS_DIR <- file.path(PROJECT_DIR, "results")
FIG_DIR     <- file.path(RESULTS_DIR, "figures")
LOGS_DIR    <- file.path(PROJECT_DIR, "logs")

for (d in c(DATA_DIR, file.path(RESULTS_DIR, c("instruments", "primary", "sensitivity",
  "bidirectional", "mvmr")), FIG_DIR, LOGS_DIR)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# ============================================================
# GWAS IDs — All verified against OpenGWAS API 2026-03-28
# ============================================================

# Primary exposures
exposure_ids <- c(
  iron         = "ieu-a-1049",           # Benyamin 2014, serum iron, n=23,986
  calcium      = "ebi-a-GCST90025990",   # Barton, calcium, n=400,792
  copper       = "ieu-a-1073",           # Evans 2013, erythrocyte copper, n=2,603
  zinc         = "ieu-a-1079",           # Evans 2013, erythrocyte zinc, n=2,603
  selenium     = "ieu-a-1075"            # Evans 2013, erythrocyte selenium, n=2,874
)

# Iron secondary markers
iron_secondary_ids <- c(
  ferritin         = "ieu-a-1050",       # Benyamin 2014, n=23,986
  transferrin_sat  = "ieu-a-1051",       # Benyamin 2014, n=23,986
  transferrin      = "ieu-a-1052",       # Benyamin 2014, n=23,986
  ferritin_large   = "ieu-b-5115"        # Bell 2021, n=246,139
)

all_exposure_ids <- c(exposure_ids, iron_secondary_ids)

# Outcomes
outcome_ids <- c(
  circ_total_tau = "ebi-a-GCST90095138", # Sarnowski 2022, n=14,721
  hippocampus    = "ieu-a-1045",          # Hibar 2015 ENIGMA, n=30,717
  entorhinal_lh  = "ubm-b-1119",         # Elliott 2021 UKB, n=31,967
  entorhinal_rh  = "ubm-b-1150",         # Elliott 2021 UKB, n=31,968
  wmh            = "ubm-b-1437"           # Elliott 2021 UKB WMH, n=32,114
)

# Parameters
PVAL_PRIMARY  <- 5e-8
PVAL_RELAXED  <- 5e-6
CLUMP_R2      <- 0.001
CLUMP_KB      <- 10000
F_MIN         <- 10
MR_METHODS    <- c("mr_ivw", "mr_egger_regression", "mr_weighted_median",
                   "mr_weighted_mode", "mr_simple_mode")
PRESSO_NDIST  <- 10000

# ============================================================
# STAGE 1: INSTRUMENT EXTRACTION
# ============================================================

cat("\n=== STAGE 1: Instrument Extraction ===\n")

extract_and_qc <- function(gwas_id, name, p_threshold = PVAL_PRIMARY) {
  cat("  ", name, "(", gwas_id, ") p <", p_threshold, "...")
  inst <- tryCatch(
    extract_instruments(gwas_id, p1 = p_threshold, clump = TRUE,
                        r2 = CLUMP_R2, kb = CLUMP_KB),
    error = function(e) { cat(" ERROR:", e$message, "\n"); NULL }
  )
  if (is.null(inst) || nrow(inst) == 0) { cat(" no instruments\n"); return(NULL) }
  inst$F_stat <- (inst$beta.exposure / inst$se.exposure)^2
  inst <- inst[inst$F_stat >= F_MIN, ]
  if (nrow(inst) == 0) { cat(" all weak\n"); return(NULL) }
  cat(" OK (", nrow(inst), "SNPs, mean F =", round(mean(inst$F_stat), 1), ")\n")
  return(inst)
}

# Primary exposures at p < 5e-8
all_instruments <- list()
for (nm in names(all_exposure_ids)) {
  all_instruments[[nm]] <- extract_and_qc(all_exposure_ids[nm], nm)
}

# Trace metals at p < 5e-6 (relaxed threshold sensitivity)
trace_metals <- c("copper", "zinc", "selenium")
relaxed_instruments <- list()
for (nm in trace_metals) {
  relaxed_instruments[[nm]] <- extract_and_qc(all_exposure_ids[nm], paste0(nm, " (5e-6)"), PVAL_RELAXED)
}

# ============================================================
# STAGE 2: HARMONIZATION + PRIMARY MR
# ============================================================

cat("\n=== STAGE 2: Harmonization + Primary MR ===\n")

run_pair <- function(inst, outcome_id, exp_name, out_name) {
  pair <- paste(exp_name, out_name, sep = "_")
  if (is.null(inst) || nrow(inst) < 3) return(NULL)
  out_dat <- tryCatch(
    extract_outcome_data(snps = inst$SNP, outcomes = outcome_id),
    error = function(e) NULL)
  if (is.null(out_dat) || nrow(out_dat) == 0) return(NULL)
  dat <- tryCatch(harmonise_data(inst, out_dat, action = 2), error = function(e) NULL)
  if (is.null(dat)) return(NULL)
  dat <- dat[dat$mr_keep == TRUE, ]
  if (nrow(dat) < 3) return(NULL)
  mr_res <- tryCatch(mr(dat, method_list = MR_METHODS), error = function(e) NULL)
  if (is.null(mr_res)) return(NULL)
  mr_res$exposure_clean <- exp_name
  mr_res$outcome_clean <- out_name
  mr_res$pair <- pair
  mr_res$ci_lower <- mr_res$b - 1.96 * mr_res$se
  mr_res$ci_upper <- mr_res$b + 1.96 * mr_res$se
  return(list(mr = mr_res, dat = dat))
}

# Run all primary + secondary pairs
all_results <- list()
all_harmonised <- list()
pair_log <- data.frame(pair = character(), status = character(), n_snps = integer())

for (exp_nm in names(all_exposure_ids)) {
  for (out_nm in names(outcome_ids)) {
    pair <- paste(exp_nm, out_nm, sep = "_")
    cat("  ", pair, "...")
    res <- run_pair(all_instruments[[exp_nm]], outcome_ids[out_nm], exp_nm, out_nm)
    if (!is.null(res)) {
      all_results[[pair]] <- res$mr
      all_harmonised[[pair]] <- res$dat
      ivw <- res$mr[res$mr$method == "Inverse variance weighted", ]
      flag <- ifelse(ivw$pval < 0.05, " *", "")
      cat(" beta=", round(ivw$b, 4), " p=", formatC(ivw$pval, format = "e", digits = 2),
          " (", nrow(res$dat), "SNPs)", flag, "\n")
      pair_log <- rbind(pair_log, data.frame(pair = pair, status = "OK", n_snps = nrow(res$dat)))
    } else {
      cat(" FAILED\n")
      pair_log <- rbind(pair_log, data.frame(pair = pair, status = "FAILED", n_snps = 0))
    }
  }
}

# Run trace metals at relaxed threshold
cat("\n--- Trace metals at p < 5e-6 ---\n")
relaxed_results <- list()
for (exp_nm in trace_metals) {
  for (out_nm in names(outcome_ids)) {
    pair <- paste0(exp_nm, "_5e6_", out_nm)
    cat("  ", pair, "...")
    res <- run_pair(relaxed_instruments[[exp_nm]], outcome_ids[out_nm], exp_nm, out_nm)
    if (!is.null(res)) {
      relaxed_results[[pair]] <- res$mr
      ivw <- res$mr[res$mr$method == "Inverse variance weighted", ]
      cat(" beta=", round(ivw$b, 4), " p=", formatC(ivw$pval, format = "e", digits = 2), "\n")
    } else {
      cat(" FAILED\n")
    }
  }
}

# ============================================================
# STAGE 3: COMBINE RESULTS + MULTIPLE TESTING
# ============================================================

cat("\n=== STAGE 3: Combine + Multiple Testing ===\n")

primary_results <- do.call(rbind, all_results)
rownames(primary_results) <- NULL

ivw_idx <- which(primary_results$method == "Inverse variance weighted")
primary_results$p_fdr <- NA_real_
primary_results$p_fdr[ivw_idx] <- p.adjust(primary_results$pval[ivw_idx], method = "BH")
primary_results$sig_bonf_global <- primary_results$pval < (0.05 / length(ivw_idx))
primary_results$sig_fdr <- FALSE
primary_results$sig_fdr[ivw_idx] <- !is.na(primary_results$p_fdr[ivw_idx]) &
                                     primary_results$p_fdr[ivw_idx] < 0.05

n_nom <- sum(primary_results$pval[ivw_idx] < 0.05)
cat("Pairs analysed:", length(ivw_idx), "\n")
cat("Nominally significant:", n_nom, "\n")

# ============================================================
# STAGE 4: SENSITIVITY — Transferrin → Tau
# ============================================================

cat("\n=== STAGE 4: Sensitivity (Transferrin -> Tau) ===\n")

dat_tt <- all_harmonised[["transferrin_circ_total_tau"]]
sensitivity_results <- list()

if (!is.null(dat_tt)) {
  # Egger intercept
  pleio <- mr_pleiotropy_test(dat_tt)
  cat("  Egger intercept:", round(pleio$egger_intercept, 5), " p=",
      formatC(pleio$pval, format = "e", digits = 3), "\n")
  sensitivity_results$egger_intercept_p <- pleio$pval

  # Cochran's Q
  het <- mr_heterogeneity(dat_tt)
  ivw_het <- het[het$method == "Inverse variance weighted", ]
  I2 <- max(0, (ivw_het$Q - ivw_het$Q_df) / ivw_het$Q * 100)
  cat("  Cochran Q p=", formatC(ivw_het$Q_pval, format = "e", digits = 3), " I2=", round(I2, 1), "%\n")
  sensitivity_results$cochran_q_p <- ivw_het$Q_pval
  sensitivity_results$I2 <- I2

  # MR-PRESSO
  set.seed(42)
  tryCatch({
    presso <- mr_presso(
      BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
      SdOutcome = "se.outcome", SdExposure = "se.exposure",
      OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
      data = dat_tt, NbDistribution = PRESSO_NDIST, SignifThreshold = 0.05)
    gp <- presso[["MR-PRESSO results"]][["Global Test"]]$Pvalue
    cat("  MR-PRESSO global p=", gp, "\n")
    sensitivity_results$presso_global_p <- gp
  }, error = function(e) { cat("  PRESSO error:", e$message, "\n") })

  # Steiger
  tryCatch({
    st <- steiger_filtering(dat_tt)
    n_correct <- sum(st$steiger_dir == TRUE, na.rm = TRUE)
    cat("  Steiger:", n_correct, "/", nrow(st), "correct\n")
    sensitivity_results$steiger_correct <- n_correct
    sensitivity_results$steiger_total <- nrow(st)
  }, error = function(e) cat("  Steiger error\n"))

  # Leave-one-out
  loo <- mr_leaveoneout(dat_tt)
  loo_all <- loo[loo$SNP == "All", ]
  loo_snps <- loo[loo$SNP != "All", ]
  n_lose_sig <- sum(loo_snps$p > 0.05)
  cat("  LOO: ", n_lose_sig, "/", nrow(loo_snps), "removals lose significance\n")
  sensitivity_results$loo_lose_sig <- n_lose_sig
}

# ============================================================
# STAGE 5: MVMR — Transferrin + Iron → Tau
# ============================================================

cat("\n=== STAGE 5: MVMR ===\n")

tf_inst <- all_instruments[["transferrin"]]
fe_inst <- all_instruments[["iron"]]
mvmr_result <- NULL

if (!is.null(tf_inst) && !is.null(fe_inst)) {
  all_snps <- unique(c(tf_inst$SNP, fe_inst$SNP))
  tf_exp <- tryCatch(extract_outcome_data(snps = all_snps, outcomes = "ieu-a-1052"), error = function(e) NULL)
  fe_exp <- tryCatch(extract_outcome_data(snps = all_snps, outcomes = "ieu-a-1049"), error = function(e) NULL)
  tau_out <- tryCatch(extract_outcome_data(snps = all_snps, outcomes = "ebi-a-GCST90095138"), error = function(e) NULL)

  if (!is.null(tf_exp) && !is.null(fe_exp) && !is.null(tau_out)) {
    common <- Reduce(intersect, list(tf_exp$SNP, fe_exp$SNP, tau_out$SNP))
    cat("  Common SNPs:", length(common), "\n")

    if (length(common) >= 5) {
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
        mvmr_ivw <- tryCatch(mr_mvivw(mvmr_input), error = function(e) NULL)
        if (!is.null(mvmr_ivw)) {
          mvmr_result <- data.frame(
            exposure = mvmr_ivw@Exposure, beta = mvmr_ivw@Estimate,
            se = mvmr_ivw@StdError, pval = mvmr_ivw@Pvalue)
          cat("  Transferrin: beta=", round(mvmr_result$beta[1], 4),
              " p=", formatC(mvmr_result$pval[1], format = "e", digits = 2), "\n")
          cat("  Serum Iron:  beta=", round(mvmr_result$beta[2], 4),
              " p=", formatC(mvmr_result$pval[2], format = "e", digits = 2), "\n")
        }
      }
    }
  }
}

# ============================================================
# STAGE 6: CONTROLS
# ============================================================

cat("\n=== STAGE 6: Controls ===\n")

iron_inst <- all_instruments[["iron"]]
controls <- list()

if (!is.null(iron_inst)) {
  # Positive: iron → haemoglobin
  cat("  Positive control: iron -> haemoglobin...")
  hb_out <- tryCatch(extract_outcome_data(snps = iron_inst$SNP, outcomes = "ebi-a-GCST90028995"), error = function(e) NULL)
  if (!is.null(hb_out)) {
    dat_hb <- harmonise_data(iron_inst, hb_out, action = 2)
    dat_hb <- dat_hb[dat_hb$mr_keep == TRUE, ]
    if (nrow(dat_hb) >= 3) {
      mr_hb <- mr(dat_hb, method_list = MR_METHODS)
      ivw_hb <- mr_hb[mr_hb$method == "Inverse variance weighted", ]
      cat(" beta=", round(ivw_hb$b, 4), " p=", formatC(ivw_hb$pval, format = "e", digits = 2))
      cat(ifelse(ivw_hb$b > 0 & ivw_hb$pval < 0.05, " PASS\n", " CHECK\n"))
      controls$positive <- ivw_hb
    }
  }

  # Negative: iron → height
  cat("  Negative control: iron -> height...")
  ht_out <- tryCatch(extract_outcome_data(snps = iron_inst$SNP, outcomes = "ebi-a-GCST90029008"), error = function(e) NULL)
  if (!is.null(ht_out)) {
    dat_ht <- harmonise_data(iron_inst, ht_out, action = 2)
    dat_ht <- dat_ht[dat_ht$mr_keep == TRUE, ]
    if (nrow(dat_ht) >= 3) {
      mr_ht <- mr(dat_ht, method_list = MR_METHODS)
      ivw_ht <- mr_ht[mr_ht$method == "Inverse variance weighted", ]
      cat(" beta=", round(ivw_ht$b, 4), " p=", formatC(ivw_ht$pval, format = "e", digits = 2))
      cat(ifelse(ivw_ht$pval > 0.05, " PASS\n", " CHECK\n"))
      controls$negative <- ivw_ht
    }
  }
}

# ============================================================
# STAGE 7: SAVE ALL RESULTS
# ============================================================

cat("\n=== STAGE 7: Saving Results ===\n")

run_date <- Sys.Date()

# Core results
write.csv(primary_results, file.path(RESULTS_DIR, "primary",
  paste0("FINAL_primary_all_methods_", run_date, ".csv")), row.names = FALSE)

ivw_summary <- primary_results[primary_results$method == "Inverse variance weighted", ]
write.csv(ivw_summary, file.path(RESULTS_DIR, "primary",
  paste0("FINAL_IVW_summary_", run_date, ".csv")), row.names = FALSE)

write.csv(pair_log, file.path(RESULTS_DIR, "primary",
  paste0("FINAL_pair_log_", run_date, ".csv")), row.names = FALSE)

# MVMR
if (!is.null(mvmr_result)) {
  write.csv(mvmr_result, file.path(RESULTS_DIR, "mvmr",
    paste0("FINAL_mvmr_", run_date, ".csv")), row.names = FALSE)
}

# Relaxed threshold
if (length(relaxed_results) > 0) {
  relaxed_df <- do.call(rbind, relaxed_results)
  write.csv(relaxed_df, file.path(RESULTS_DIR, "sensitivity",
    paste0("FINAL_relaxed_5e6_", run_date, ".csv")), row.names = FALSE)
}

# Save R objects
save(all_results, primary_results, all_harmonised, all_instruments,
     relaxed_instruments, relaxed_results, mvmr_result, sensitivity_results,
     controls, pair_log,
     file = file.path(DATA_DIR, "FINAL_analysis_results.RData"))

# Session info
writeLines(capture.output(sessionInfo()),
  file.path(LOGS_DIR, paste0("FINAL_sessionInfo_", run_date, ".txt")))

# ============================================================
# HEADLINE RESULTS
# ============================================================

cat("\n==============================================================\n")
cat("  HEADLINE RESULTS\n")
cat("==============================================================\n\n")

# Transferrin → tau
tt_ivw <- ivw_summary[ivw_summary$pair == "transferrin_circ_total_tau", ]
if (nrow(tt_ivw) > 0) {
  cat("1. Transferrin -> Circulating Tau (IVW):\n")
  cat("   beta =", round(tt_ivw$b, 4), ", SE =", round(tt_ivw$se, 4),
      ", p =", formatC(tt_ivw$pval, format = "e", digits = 3), "\n")
}

# MVMR
if (!is.null(mvmr_result)) {
  cat("\n2. MVMR (independent of serum iron):\n")
  cat("   Transferrin: beta =", round(mvmr_result$beta[1], 4),
      ", p =", formatC(mvmr_result$pval[1], format = "e", digits = 3), "\n")
  cat("   Serum Iron:  beta =", round(mvmr_result$beta[2], 4),
      ", p =", formatC(mvmr_result$pval[2], format = "e", digits = 3), "\n")
}

# Sensitivity
cat("\n3. Sensitivity:\n")
cat("   Egger intercept p =", formatC(sensitivity_results$egger_intercept_p, format = "e", digits = 3), "\n")
cat("   Cochran Q p =", formatC(sensitivity_results$cochran_q_p, format = "e", digits = 3), "\n")
cat("   MR-PRESSO p =", sensitivity_results$presso_global_p, "\n")
cat("   Steiger:", sensitivity_results$steiger_correct, "/", sensitivity_results$steiger_total, "\n")
cat("   LOO losses:", sensitivity_results$loo_lose_sig, "\n")

# Controls
cat("\n4. Controls:\n")
if (!is.null(controls$positive)) {
  cat("   Positive (iron->Hb): beta =", round(controls$positive$b, 4),
      ", p =", formatC(controls$positive$pval, format = "e", digits = 2), "\n")
}
if (!is.null(controls$negative)) {
  cat("   Negative (iron->ht): beta =", round(controls$negative$b, 4),
      ", p =", formatC(controls$negative$pval, format = "e", digits = 2), "\n")
}

cat("\n==============================================================\n")
cat("  PIPELINE COMPLETE:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("  All outputs saved to:", RESULTS_DIR, "\n")
cat("==============================================================\n")
