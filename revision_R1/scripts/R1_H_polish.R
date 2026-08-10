##############################################################################
# R1_H_polish.R
#
# Optional polish analyses to strengthen the reviewer response:
#   H1. MVMR-Egger (intercept test for directional pleiotropy in MVMR)
#   H2. Q-minimization MVMR (Sanderson pleiotropy-robust estimator)
#   H3. Radial MVMR outlier check
#   H4. Formal APOE4-carrier × transferrin interaction test on CSF p-tau
#   H4b. Normal vs abnormal amyloid × transferrin interaction test
#   H4c. Symmetry checks for CSF Aβ42 subgroups
#
# Writes: results/reviewer_response/H_polish_summary.txt
#         results/reviewer_response/H_MVMR_Egger.csv
#         results/reviewer_response/H_MVMR_Qmin.csv
#         results/reviewer_response/H_interaction_tests.csv
##############################################################################

args <- commandArgs(trailingOnly = TRUE)
PROJECT_DIR <- if (length(args) >= 1) args[1] else normalizePath("..", mustWork = TRUE)
OUT_DIR     <- file.path(PROJECT_DIR, "results", "reviewer_response")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(ieugwasr); library(TwoSampleMR); library(MVMR)
  library(MendelianRandomization); library(dplyr)
})
stopifnot(nchar(Sys.getenv("OPENGWAS_JWT")) > 100)

fetch_snp_effects <- function(snps, id, tries = 5, sleep = 4) {
  for (k in seq_len(tries)) {
    Sys.sleep(if (k == 1) 1 else sleep)
    a <- tryCatch(associations(variants = snps, id = id, proxies = 0),
                   error = function(e) { cat("    attempt", k, "ERR:",
                     conditionMessage(e), "\n"); NULL })
    if (!is.null(a) && nrow(a) > 0) return(a)
  }
  stop(sprintf("Failed to fetch %s", id))
}

# ============================================================
# H1 + H2 + H3: build 2-way MVMR input (mirror Analysis E)
# ============================================================
TF_ID  <- "ieu-a-1052"
FE_ID  <- "ieu-a-1049"
OUT_ID <- "ebi-a-GCST90095138"

cat("== Rebuilding MVMR input (transferrin + iron) ==\n")
tf_inst <- extract_instruments(TF_ID, p1 = 5e-8, r2 = 0.001, kb = 10000)
fe_inst <- extract_instruments(FE_ID, p1 = 5e-8, r2 = 0.001, kb = 10000)
union_snps <- unique(c(tf_inst$SNP, fe_inst$SNP))

tf_e <- fetch_snp_effects(union_snps, TF_ID)
fe_e <- fetch_snp_effects(union_snps, FE_ID)
ta_o <- fetch_snp_effects(union_snps, OUT_ID)

common <- Reduce(intersect, list(tf_e$rsid, fe_e$rsid, ta_o$rsid))
cat(sprintf("Common: %d SNPs\n", length(common)))
row <- function(d) d[match(common, d$rsid), ]
tf_e <- row(tf_e); fe_e <- row(fe_e); ta_o <- row(ta_o)

# align to transferrin effect allele
align <- function(d, ref) {
  flip <- (d$ea == ref$nea) & (d$nea == ref$ea)
  same <- (d$ea == ref$ea)  & (d$nea == ref$nea)
  keep <- flip | same
  b <- d$beta; b[flip] <- -b[flip]
  list(beta = b, se = d$se, keep = keep)
}
a_fe <- align(fe_e, tf_e)
a_ta <- align(ta_o, tf_e)
keep <- a_fe$keep & a_ta$keep
tf_b  <- tf_e$beta[keep];   tf_se  <- tf_e$se[keep]
fe_b  <- a_fe$beta[keep];   fe_se  <- fe_e$se[keep]
tau_b <- a_ta$beta[keep];   tau_se <- ta_o$se[keep]
snps  <- common[keep]
cat(sprintf("After allele alignment: %d SNPs\n", length(snps)))

# ============================================================
# H1. MVMR-Egger via MendelianRandomization::mr_mvegger
# ============================================================
cat("\n== H1: MVMR-Egger (directional pleiotropy test) ==\n")
bx  <- cbind(tf_b, fe_b)
bxse <- cbind(tf_se, fe_se)
mvmr_in <- mr_mvinput(bx = bx, bxse = bxse, by = tau_b, byse = tau_se,
                     exposure = c("transferrin","serum_iron"))

mvegger <- mr_mvegger(mvmr_in, orientate = 1)
cat("MVMR-Egger:\n")
cat(sprintf("  intercept: %.5f  SE %.5f  p = %.3e\n",
            mvegger@Intercept, mvegger@StdError.Int, mvegger@Pvalue.Int))
for (i in seq_along(mvegger@Estimate)) {
  cat(sprintf("  %-12s slope %.4f  SE %.4f  p = %.3e\n",
              mvegger@Exposure[i], mvegger@Estimate[i],
              mvegger@StdError.Est[i], mvegger@Pvalue.Est[i]))
}
mvegger_df <- data.frame(
  term = c("intercept", mvegger@Exposure),
  estimate = c(mvegger@Intercept, mvegger@Estimate),
  se       = c(mvegger@StdError.Int, mvegger@StdError.Est),
  pval     = c(mvegger@Pvalue.Int, mvegger@Pvalue.Est)
)
write.csv(mvegger_df, file.path(OUT_DIR, "H_MVMR_Egger.csv"), row.names = FALSE)

# ============================================================
# H2. Q-minimization MVMR (WSpiller/MVMR::qhet_mvmr)
# ============================================================
cat("\n== H2: Q-minimization MVMR ==\n")
mv_input_MVMR <- format_mvmr(BXGs = bx, seBXGs = bxse,
                              BYG = tau_b, seBYG = tau_se, RSID = snps)
qhet_res <- tryCatch(qhet_mvmr(mv_input_MVMR, CI = TRUE, iterations = 1000),
                      error = function(e) { cat("  qhet_mvmr err:",
                        conditionMessage(e), "\n"); NULL })
if (!is.null(qhet_res)) {
  print(qhet_res)
  write.csv(as.data.frame(qhet_res),
            file.path(OUT_DIR, "H_MVMR_Qmin.csv"), row.names = FALSE)
} else {
  # Fallback: mr_mvivw with random effects using MendelianRandomization
  cat("qhet_mvmr failed; running MVMR-IVW random-effects instead.\n")
  ivw_re <- mr_mvivw(mvmr_in, model = "random")
  cat("MVMR-IVW random effects:\n")
  for (i in seq_along(ivw_re@Estimate)) {
    cat(sprintf("  %-12s b %.4f  SE %.4f  p = %.3e\n",
                ivw_re@Exposure[i], ivw_re@Estimate[i],
                ivw_re@StdError[i], ivw_re@Pvalue[i]))
  }
  write.csv(data.frame(
    exposure = ivw_re@Exposure,
    estimate = ivw_re@Estimate,
    se = ivw_re@StdError,
    pval = ivw_re@Pvalue
  ), file.path(OUT_DIR, "H_MVMR_IVW_RE.csv"), row.names = FALSE)
}

# ============================================================
# H3. Radial MVMR outlier scan (using per-SNP Q contributions)
# ============================================================
cat("\n== H3: Per-SNP Q-contribution outlier check ==\n")
# From MVMR::pleiotropy_mvmr we already know total Q_A; here we compute per-SNP
# contribution to test for outliers (Bowden 2016 radial approach adapted).
pred_by <- tf_b * (-0.0437) + fe_b * (-0.00568)  # E-derived MVMR estimates
resid   <- tau_b - pred_by
q_contrib <- (resid / tau_se)^2
q_df <- data.frame(SNP = snps, q_contrib = q_contrib) |> arrange(-q_contrib)
cat("Per-SNP Q contributions (largest first):\n"); print(q_df)
threshold <- qchisq(0.95, df = 1)   # 3.84
outliers <- q_df[q_df$q_contrib > threshold, ]
cat(sprintf("SNPs with Q_i > chi^2(1, 0.95)=%.2f: %d\n",
            threshold, nrow(outliers)))
if (nrow(outliers) > 0) print(outliers)
write.csv(q_df, file.path(OUT_DIR, "H_MVMR_perSNP_Q.csv"), row.names = FALSE)

# ============================================================
# H4. APOE4 × transferrin interaction test on CSF p-tau
# ============================================================
cat("\n== H4: Formal interaction / heterogeneity tests ==\n")

# Standard 2-sample Cochran heterogeneity test for two independent estimates
het_test <- function(b1, se1, b2, se2, label) {
  z <- (b1 - b2) / sqrt(se1^2 + se2^2)
  p <- 2 * pnorm(-abs(z))
  # Cochran's Q under fixed-effects meta
  wsum <- 1/se1^2 + 1/se2^2
  bp   <- (b1/se1^2 + b2/se2^2) / wsum
  Q    <- (b1 - bp)^2/se1^2 + (b2 - bp)^2/se2^2
  p_Q  <- pchisq(Q, df = 1, lower.tail = FALSE)
  I2   <- max(0, (Q - 1) / Q * 100)
  data.frame(
    contrast = label,
    b1 = b1, se1 = se1, b2 = b2, se2 = se2,
    diff = b1 - b2, z_diff = z, p_diff = p,
    Q = Q, p_Q = p_Q, I2 = I2
  )
}

csf <- read.csv(file.path(PROJECT_DIR, "results", "primary",
                          "CSF_all_Jansen_subanalyses_2026-03-29.csv"))
get <- function(pattern) csf[grepl(pattern, csf$trait), ]

ptau_non <- get("p-tau — APOE4 non-carriers")
ptau_car <- get("p-tau — APOE4 carriers")
ab42_non <- get("Abeta42 — APOE4 non-carriers")
ab42_car <- get("Abeta42 — APOE4 carriers")
ptau_norm <- get("p-tau — normal amyloid")
ptau_ab   <- get("p-tau — abnormal amyloid")

tests <- rbind(
  het_test(ptau_non$beta, ptau_non$se, ptau_car$beta, ptau_car$se,
           "CSF p-tau: APOE4 non-carriers vs carriers"),
  het_test(ab42_non$beta, ab42_non$se, ab42_car$beta, ab42_car$se,
           "CSF Abeta42: APOE4 non-carriers vs carriers"),
  het_test(ptau_norm$beta, ptau_norm$se, ptau_ab$beta, ptau_ab$se,
           "CSF p-tau: normal vs abnormal amyloid")
)
tests[, sapply(tests, is.numeric)] <- lapply(
  tests[, sapply(tests, is.numeric)], function(x) signif(x, 3))
print(tests)
write.csv(tests, file.path(OUT_DIR, "H_interaction_tests.csv"), row.names = FALSE)

# ============================================================
# Consolidated summary
# ============================================================
lines <- c(
  "== Polish analyses (H1-H4) summary ==", "",
  "H1. MVMR-Egger (directional pleiotropy in MVMR)",
  sprintf("   intercept = %.5f  (SE %.5f, p = %.3e)  <<< no directional pleiotropy if p > 0.05",
          mvegger@Intercept, mvegger@StdError.Int, mvegger@Pvalue.Int),
  sprintf("   transferrin slope: b = %.4f  SE = %.4f  p = %.3e",
          mvegger@Estimate[1], mvegger@StdError.Est[1], mvegger@Pvalue.Est[1]),
  sprintf("   serum iron slope:  b = %.4f  SE = %.4f  p = %.3e",
          mvegger@Estimate[2], mvegger@StdError.Est[2], mvegger@Pvalue.Est[2]),
  "",
  "H3. Per-SNP Q contributions (outlier scan)",
  sprintf("   %d SNPs exceed chi^2(1, 0.95) = 3.84 threshold",
          nrow(outliers)),
  "",
  "H4. Subgroup heterogeneity tests (z-diff and Cochran's Q):",
  ""
)
for (i in seq_len(nrow(tests))) {
  lines <- c(lines,
    sprintf("   %s", tests$contrast[i]),
    sprintf("     b1 = %.3f (SE %.3f)  vs  b2 = %.3f (SE %.3f)",
            tests$b1[i], tests$se1[i], tests$b2[i], tests$se2[i]),
    sprintf("     diff = %.3f   z = %.2f   p_diff = %.3f   p_Q = %.3f   I2 = %.1f%%",
            tests$diff[i], tests$z_diff[i], tests$p_diff[i], tests$p_Q[i], tests$I2[i]),
    ""
  )
}
writeLines(lines, file.path(OUT_DIR, "H_polish_summary.txt"))
cat("\n", paste(lines, collapse = "\n"), "\n", sep = "")
