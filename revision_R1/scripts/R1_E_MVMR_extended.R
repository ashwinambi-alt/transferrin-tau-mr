##############################################################################
# R1_E_MVMR_extended.R
#
# Reviewer response E: Extended MVMR reporting.
#
# Re-runs 2-exposure MVMR (transferrin + serum iron -> circulating total-tau)
# and reports what the paper currently omits:
#   - explicit instrument union set
#   - Sanderson-Windmeijer conditional F-statistics for each exposure
#   - Q_A heterogeneity
#   - exposure-exposure correlation among the instruments
#   - overlap disclosure (Benyamin transferrin and iron share same cohort)
#
# Uses the MVMR package (WSpiller/MVMR).
##############################################################################

args <- commandArgs(trailingOnly = TRUE)
PROJECT_DIR <- if (length(args) >= 1) args[1] else normalizePath("..", mustWork = TRUE)
OUT_DIR     <- file.path(PROJECT_DIR, "results", "reviewer_response")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(ieugwasr); library(TwoSampleMR); library(MVMR); library(dplyr)
})
stopifnot(nchar(Sys.getenv("OPENGWAS_JWT")) > 100)

# TwoSampleMR extract_outcome_data was hitting 401 despite valid JWT — switch
# to ieugwasr::associations directly and wrap into TwoSampleMR format ourselves.
fetch_snp_effects <- function(snps, id, tries = 5, sleep = 4) {
  for (k in seq_len(tries)) {
    Sys.sleep(if (k == 1) 1 else sleep)
    a <- tryCatch(
      associations(variants = snps, id = id, proxies = 0),
      error = function(e) { cat(sprintf("    attempt %d ERR: %s\n", k, conditionMessage(e))); NULL }
    )
    if (!is.null(a) && nrow(a) > 0) {
      return(data.frame(
        SNP = a$rsid,
        beta.outcome = a$beta,
        se.outcome   = a$se,
        effect_allele.outcome = a$ea,
        other_allele.outcome  = a$nea,
        eaf.outcome  = if ("eaf" %in% names(a)) a$eaf else NA,
        pval.outcome = a$p,
        samplesize.outcome = if ("n" %in% names(a)) a$n else NA,
        outcome = id, id.outcome = id, mr_keep = TRUE,
        stringsAsFactors = FALSE
      ))
    }
  }
  stop(sprintf("Failed to fetch %s after %d attempts", id, tries))
}

TF_ID  <- "ieu-a-1052"          # transferrin, Benyamin 2014 (n=23,986)
FE_ID  <- "ieu-a-1049"          # serum iron,  Benyamin 2014 (n=23,986) — SAME cohort
OUT_ID <- "ebi-a-GCST90095138"  # total tau,   Sarnowski 2022 (n=14,721)

# ---- 1. Build union instrument set ----
cat("== Extended MVMR: transferrin + serum iron -> total-tau ==\n")
cat("Step 1: pull genome-wide-significant instruments for each exposure\n")

tf_inst <- extract_instruments(TF_ID, p1 = 5e-8, r2 = 0.001, kb = 10000)
fe_inst <- extract_instruments(FE_ID, p1 = 5e-8, r2 = 0.001, kb = 10000)
cat(sprintf("  Transferrin: %d SNPs  |  Serum iron: %d SNPs\n",
            nrow(tf_inst), nrow(fe_inst)))
union_snps <- unique(c(tf_inst$SNP, fe_inst$SNP))
cat(sprintf("  Union: %d SNPs\n", length(union_snps)))

# ---- 2. Extract each exposure's SNP-effect on the union set ----
cat("Step 2: fetch exposure SNP-effects on the union set\n")
tf_exp  <- fetch_snp_effects(union_snps, TF_ID)
fe_exp  <- fetch_snp_effects(union_snps, FE_ID)
tau_out <- fetch_snp_effects(union_snps, OUT_ID)

# Intersect
common <- Reduce(intersect, list(tf_exp$SNP, fe_exp$SNP, tau_out$SNP))
cat(sprintf("  Common (transferrin + iron + tau): %d SNPs\n", length(common)))

get_row <- function(df, snp) df[match(snp, df$SNP), ]

tf_r  <- get_row(tf_exp,  common)
fe_r  <- get_row(fe_exp,  common)
tau_r <- get_row(tau_out, common)

# Align effect alleles across exposures + outcome (align to transferrin's effect allele)
align_beta <- function(target_ea, target_nea, target_beta,
                       ref_ea,    ref_nea) {
  # If ref effect allele matches target other allele, flip sign of target beta
  flip <- (target_ea == ref_nea) & (target_nea == ref_ea)
  same <- (target_ea == ref_ea)  & (target_nea == ref_nea)
  b <- target_beta
  b[flip] <- -b[flip]
  keep <- flip | same
  list(beta = b, keep = keep)
}

# align fe and tau to tf
a_fe  <- align_beta(fe_r$effect_allele.outcome,  fe_r$other_allele.outcome,
                    fe_r$beta.outcome,
                    tf_r$effect_allele.outcome,  tf_r$other_allele.outcome)
a_tau <- align_beta(tau_r$effect_allele.outcome, tau_r$other_allele.outcome,
                    tau_r$beta.outcome,
                    tf_r$effect_allele.outcome,  tf_r$other_allele.outcome)

keep <- a_fe$keep & a_tau$keep
cat(sprintf("  Allele-alignable SNPs: %d / %d\n", sum(keep), length(common)))
tf_r  <- tf_r[keep, ];  fe_r  <- fe_r[keep, ]; tau_r <- tau_r[keep, ]
a_fe$beta  <- a_fe$beta[keep]
a_tau$beta <- a_tau$beta[keep]
common <- common[keep]

# ---- 3. Build MVMR input ----
bx  <- cbind(tf_r$beta.outcome, a_fe$beta)
bxse <- cbind(tf_r$se.outcome,  fe_r$se.outcome)
by   <- a_tau$beta
byse <- tau_r$se.outcome
colnames(bx) <- colnames(bxse) <- c("transferrin", "serum_iron")

mv_input <- format_mvmr(BXGs = bx, seBXGs = bxse, BYG = by, seBYG = byse,
                        RSID = common)

# ---- 4. Sanderson-Windmeijer conditional F-statistics ----
# Uses MVMR::strength_mvmr assuming no covariance (independent exposure GWAS estimates
# from summary statistics; conservative given Benyamin cohort overlap).
Fcond <- strength_mvmr(r_input = mv_input, gencov = 0)
cat("\nSanderson-Windmeijer conditional F-statistics:\n"); print(Fcond)

# ---- 5. MVMR IVW + Q_A heterogeneity ----
mv_res <- ivw_mvmr(r_input = mv_input)
cat("\nMVMR IVW:\n"); print(mv_res)

Q_A <- pleiotropy_mvmr(r_input = mv_input, gencov = 0)
cat("\nMVMR Q_A heterogeneity:\n"); print(Q_A)

# ---- 6. Exposure correlation ----
cat("\nPearson correlation between transferrin and serum iron SNP-effects (aligned):\n")
r_exp <- cor(bx[, "transferrin"], bx[, "serum_iron"])
cat(sprintf("  r = %.3f (n = %d SNPs)\n", r_exp, nrow(bx)))

# ---- 7. Save ----
write.csv(data.frame(SNP = common,
                      beta_tf = bx[,1], se_tf = bxse[,1],
                      beta_fe = bx[,2], se_fe = bxse[,2],
                      beta_tau = by,    se_tau = byse),
          file.path(OUT_DIR, "E_MVMR_extended_snps.csv"), row.names = FALSE)

mvres_df <- data.frame(
  exposure = c("transferrin", "serum_iron"),
  beta     = mv_res[, "Estimate"],
  se       = mv_res[, "Std. Error"],
  z        = mv_res[, "t value"],
  pval     = mv_res[, "Pr(>|t|)"]
)
# MVMR::strength_mvmr labels columns "exposure1" ... "exposureK"; use position
mvres_df$cond_F <- as.numeric(Fcond[1, ])
write.csv(mvres_df, file.path(OUT_DIR, "E_MVMR_extended.csv"), row.names = FALSE)

# ---- 8. Human summary ----
lines <- c(
  "== Analysis E: Extended MVMR (transferrin + serum iron -> total-tau) ==",
  sprintf("Union SNP set (pre-alignment): %d", length(unique(c(tf_inst$SNP, fe_inst$SNP)))),
  sprintf("Retained after allele alignment and intersection with outcome: %d", nrow(bx)),
  "",
  "MVMR IVW estimates:",
  sprintf("  Transferrin: b = %.4f, SE = %.4f, p = %.3e",
          mvres_df$beta[1], mvres_df$se[1], mvres_df$pval[1]),
  sprintf("  Serum iron:  b = %.4f, SE = %.4f, p = %.3e",
          mvres_df$beta[2], mvres_df$se[2], mvres_df$pval[2]),
  "",
  "Sanderson-Windmeijer conditional F-statistics:",
  sprintf("  Transferrin cond-F: %.2f  (>10 = strong)", as.numeric(Fcond[1, 1])),
  sprintf("  Serum iron  cond-F: %.2f  (>10 = strong)", as.numeric(Fcond[1, 2])),
  "",
  sprintf("MVMR Q_A heterogeneity: Q = %.2f (p = %.3e)",
          as.numeric(Q_A$Qstat), as.numeric(Q_A$Qpval)),
  "",
  sprintf("SNP-effect Pearson correlation transferrin vs. iron: r = %.3f", r_exp),
  "",
  "Sample overlap disclosure: transferrin and serum iron come from the same",
  "Benyamin et al. 2014 cohort (n = 23,986). Total-tau (Sarnowski 2022, n = 14,721)",
  "does not overlap the Benyamin cohort. MVMR under complete exposure-exposure",
  "sample overlap is unbiased when the outcome GWAS is independent, as here."
)
writeLines(lines, file.path(OUT_DIR, "E_MVMR_extended_summary.txt"))
cat("\n", paste(lines, collapse="\n"), "\n", sep="")
