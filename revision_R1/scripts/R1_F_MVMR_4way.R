##############################################################################
# R1_F_MVMR_4way.R
#
# Reviewer response F: MVMR sensitivity with transferrin saturation and
# ferritin as additional exposures.
#
# Model: transferrin + serum_iron + ferritin + tsat -> circulating total-tau
#
# The reviewer asked to check whether transferrin's effect survives when
# ferritin (a marker of body iron STORES) and TSAT (a marker of iron
# distribution) are added alongside serum iron.
##############################################################################

args <- commandArgs(trailingOnly = TRUE)
PROJECT_DIR <- if (length(args) >= 1) args[1] else normalizePath("..", mustWork = TRUE)
OUT_DIR     <- file.path(PROJECT_DIR, "results", "reviewer_response")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(ieugwasr); library(TwoSampleMR); library(MVMR); library(dplyr)
})
stopifnot(nchar(Sys.getenv("OPENGWAS_JWT")) > 100)

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

EXP <- list(
  transferrin = "ieu-a-1052",   # Benyamin 2014
  serum_iron  = "ieu-a-1049",   # Benyamin 2014 (same cohort as transferrin)
  ferritin    = "ieu-b-5115",   # Bell 2021 (n=246,139)
  tsat        = "ieu-a-1051"    # Benyamin 2014 (same cohort)
)
OUT_ID <- "ebi-a-GCST90095138"  # Sarnowski 2022 total tau

# --- 1. Pull instruments for each exposure ---
cat("== 4-way MVMR: transferrin + iron + ferritin + TSAT -> total tau ==\n")
insts <- lapply(EXP, function(id) {
  cat(sprintf("  fetching %s...\n", id))
  extract_instruments(id, p1 = 5e-8, r2 = 0.001, kb = 10000)
})
sapply(insts, nrow) |> print()

union_snps <- unique(unlist(lapply(insts, function(x) x$SNP)))
cat(sprintf("Union: %d SNPs\n", length(union_snps)))

# --- 2. Extract each exposure's SNP-effects on the union set ---
cat("Fetching exposure and outcome effect estimates on the union set...\n")
exp_data <- lapply(EXP, function(id) fetch_snp_effects(union_snps, id))
tau_out  <- fetch_snp_effects(union_snps, OUT_ID)

common <- Reduce(intersect, c(lapply(exp_data, function(x) x$SNP),
                              list(tau_out$SNP)))
cat(sprintf("Common across all 4 exposures + outcome: %d SNPs\n", length(common)))

# Align every exposure and outcome to the transferrin effect allele
ref <- exp_data$transferrin
ref <- ref[match(common, ref$SNP), ]

align <- function(d) {
  d <- d[match(common, d$SNP), ]
  flip <- (d$effect_allele.outcome == ref$other_allele.outcome) &
          (d$other_allele.outcome  == ref$effect_allele.outcome)
  same <- (d$effect_allele.outcome == ref$effect_allele.outcome) &
          (d$other_allele.outcome  == ref$other_allele.outcome)
  b <- d$beta.outcome; b[flip] <- -b[flip]
  list(beta = b, se = d$se.outcome, keep = flip | same)
}
al <- lapply(exp_data, align)
al_tau <- align(tau_out)

keep <- Reduce(`&`, c(lapply(al, function(x) x$keep), list(al_tau$keep)))
cat(sprintf("Allele-alignable SNPs: %d / %d\n", sum(keep), length(common)))

common <- common[keep]
bx  <- do.call(cbind, lapply(al, function(x) x$beta[keep]))
bxs <- do.call(cbind, lapply(al, function(x) x$se[keep]))
by  <- al_tau$beta[keep]
bys <- al_tau$se[keep]
colnames(bx) <- colnames(bxs) <- names(EXP)

# --- 3. MVMR ---
mv_input <- format_mvmr(BXGs = bx, seBXGs = bxs, BYG = by, seBYG = bys, RSID = common)
Fcond <- strength_mvmr(r_input = mv_input, gencov = 0)
cat("\nConditional F-stats:\n"); print(Fcond)

mv_res <- ivw_mvmr(r_input = mv_input)
cat("\nMVMR IVW:\n"); print(mv_res)

Q_A <- pleiotropy_mvmr(r_input = mv_input, gencov = 0)
cat("\nQ_A:\n"); print(Q_A)

# --- 4. Also run a 3-way that drops serum_iron (since it's redundant with ferritin+TSAT) ---
cat("\n== Sensitivity: 3-way (transferrin + ferritin + TSAT) ==\n")
bx3  <- bx[, c("transferrin","ferritin","tsat")]
bxs3 <- bxs[, c("transferrin","ferritin","tsat")]
mv_input3 <- format_mvmr(BXGs = bx3, seBXGs = bxs3, BYG = by, seBYG = bys, RSID = common)
Fcond3 <- strength_mvmr(r_input = mv_input3, gencov = 0)
mv_res3 <- ivw_mvmr(r_input = mv_input3)
cat("\n3-way conditional F:\n"); print(Fcond3)
cat("\n3-way MVMR IVW:\n"); print(mv_res3)

# --- 5. Save ---
mvres_df <- data.frame(
  model    = "4-way (transferrin + iron + ferritin + tsat)",
  exposure = rownames(mv_res),
  beta     = mv_res[, "Estimate"],
  se       = mv_res[, "Std. Error"],
  pval     = mv_res[, "Pr(>|t|)"],
  cond_F   = as.numeric(Fcond[1, rownames(mv_res)])
)
mvres3_df <- data.frame(
  model    = "3-way (transferrin + ferritin + tsat)",
  exposure = rownames(mv_res3),
  beta     = mv_res3[, "Estimate"],
  se       = mv_res3[, "Std. Error"],
  pval     = mv_res3[, "Pr(>|t|)"],
  cond_F   = as.numeric(Fcond3[1, rownames(mv_res3)])
)
write.csv(rbind(mvres_df, mvres3_df),
          file.path(OUT_DIR, "F_MVMR_4way.csv"), row.names = FALSE)

lines <- c(
  "== Analysis F: 4-way MVMR sensitivity ==",
  sprintf("Instruments used: %d SNPs (allele-aligned across 4 exposures + outcome)", nrow(bx)),
  "",
  "-- 4-way (transferrin + serum iron + ferritin + TSAT) --",
  sprintf("  transferrin: b = %.4f, SE = %.4f, p = %.3e, cond-F = %.1f",
          mvres_df$beta[1], mvres_df$se[1], mvres_df$pval[1], mvres_df$cond_F[1]),
  sprintf("  serum iron : b = %.4f, SE = %.4f, p = %.3e, cond-F = %.1f",
          mvres_df$beta[2], mvres_df$se[2], mvres_df$pval[2], mvres_df$cond_F[2]),
  sprintf("  ferritin   : b = %.4f, SE = %.4f, p = %.3e, cond-F = %.1f",
          mvres_df$beta[3], mvres_df$se[3], mvres_df$pval[3], mvres_df$cond_F[3]),
  sprintf("  TSAT       : b = %.4f, SE = %.4f, p = %.3e, cond-F = %.1f",
          mvres_df$beta[4], mvres_df$se[4], mvres_df$pval[4], mvres_df$cond_F[4]),
  "",
  "-- 3-way (transferrin + ferritin + TSAT) --",
  sprintf("  transferrin: b = %.4f, SE = %.4f, p = %.3e, cond-F = %.1f",
          mvres3_df$beta[1], mvres3_df$se[1], mvres3_df$pval[1], mvres3_df$cond_F[1]),
  sprintf("  ferritin   : b = %.4f, SE = %.4f, p = %.3e, cond-F = %.1f",
          mvres3_df$beta[2], mvres3_df$se[2], mvres3_df$pval[2], mvres3_df$cond_F[2]),
  sprintf("  TSAT       : b = %.4f, SE = %.4f, p = %.3e, cond-F = %.1f",
          mvres3_df$beta[3], mvres3_df$se[3], mvres3_df$pval[3], mvres3_df$cond_F[3])
)
writeLines(lines, file.path(OUT_DIR, "F_MVMR_4way_summary.txt"))
cat("\n", paste(lines, collapse="\n"), "\n", sep="")
