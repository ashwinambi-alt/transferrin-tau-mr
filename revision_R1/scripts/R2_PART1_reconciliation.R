##############################################################################
# R2_PART1_reconciliation.R  (JAD revision round 2, Part 1)
#
# Reconciles the five numeric conflicts between manuscript text and figures.
# Base R only; no API. Reads cached harmonized data + reviewer-response CSVs.
#
# Writes:
#   results/reviewer_response/reconciliation.csv
#   results/reviewer_response/reconciliation_notes.md
##############################################################################

PROJECT_DIR <- "C:/Users/ashwi/OneDrive/Documents/EB1A Docs/MR Paper"
REV <- file.path(PROJECT_DIR, "results", "reviewer_response")
DATE_RUN <- "2026-08-07"

harm <- read.csv(file.path(PROJECT_DIR, "results", "primary",
                 "harmonised_transferrin_tau_2026-03-28.csv"),
                 stringsAsFactors = FALSE)
harm <- harm[harm$mr_keep == TRUE, ]
harm$wald    <- harm$beta.outcome / harm$beta.exposure
harm$wald_se <- harm$se.outcome / abs(harm$beta.exposure)

# ---- IVW estimators: fixed-effect (FE) and multiplicative random-effects (RE)
ivw <- function(b, s) {
  w  <- 1 / s^2
  bh <- sum(w * b) / sum(w)
  se_fe <- sqrt(1 / sum(w))
  k  <- length(b)
  Q  <- sum(w * (b - bh)^2)
  Qp <- pchisq(Q, df = k - 1, lower.tail = FALSE)
  I2 <- max(0, (Q - (k - 1)) / Q) * 100
  # multiplicative RE: inflate SE by sqrt(Q/(k-1)) but never below FE
  phi <- if (k > 1) max(1, Q / (k - 1)) else 1
  se_re <- se_fe * sqrt(phi)
  list(beta = bh, se_fe = se_fe, se_re = se_re,
       p_fe = 2 * pnorm(-abs(bh / se_fe)),
       p_re = 2 * pnorm(-abs(bh / se_re)),
       Q = Q, Q_p = Qp, I2 = I2, k = k)
}

locus <- list(TF=c("rs8177240","rs7646473","rs9990333"), HFE="rs1800562",
              HLA="rs9268633", FADS="rs174577", NAT2="rs1495741", TFRC="rs744653")

subsets <- list(
  "all_8"                          = harm$SNP,
  "excl_rs8177240"                 = setdiff(harm$SNP, "rs8177240"),
  "excl_HLA"                       = setdiff(harm$SNP, locus$HLA),
  "excl_FADS"                      = setdiff(harm$SNP, locus$FADS),
  "excl_HFE"                       = setdiff(harm$SNP, locus$HFE),
  "excl_NAT2"                      = setdiff(harm$SNP, locus$NAT2),
  "excl_HLA_FADS"                  = setdiff(harm$SNP, c(locus$HLA, locus$FADS)),
  "excl_HLA_FADS_HFE_NAT2 (joint)" = setdiff(harm$SNP, c(locus$HLA,locus$FADS,locus$HFE,locus$NAT2)),
  "TF_region_subset_3snps"         = c("rs8177240","rs9990333","rs7646473"),
  "rs8177240_alone (Wald)"         = "rs8177240"
)

rows <- list()
add <- function(analysis, n, snps, beta, se, p, method, script) {
  rows[[length(rows)+1]] <<- data.frame(
    analysis=analysis, n_snps=n, snps=snps,
    beta=signif(beta,4), se=signif(se,4), pval=signif(p,4),
    method=method, script=script, date_run=DATE_RUN, stringsAsFactors=FALSE)
}

# ---- leave-locus-out in BOTH FE and RE (resolves conflict 1.2) ----
for (nm in names(subsets)) {
  d <- harm[harm$SNP %in% subsets[[nm]], ]
  if (nrow(d) == 0) next
  snpstr <- paste(sort(d$SNP), collapse=";")
  if (nrow(d) == 1) {
    add(nm, 1, snpstr, d$wald, d$wald_se, 2*pnorm(-abs(d$wald/d$wald_se)),
        "Wald ratio", "R2_PART1_reconciliation.R")
  } else {
    r <- ivw(d$wald, d$wald_se)
    add(nm, r$k, snpstr, r$beta, r$se_fe, r$p_fe, "IVW fixed-effect",  "R2_PART1_reconciliation.R")
    add(nm, r$k, snpstr, r$beta, r$se_re, r$p_re, "IVW random-effect", "R2_PART1_reconciliation.R")
  }
}

# ---- cis-only (Analysis C) from cached files (Fig S9) ----
cis_mr  <- read.csv(file.path(REV,"C_cis_TF_mr_500kb.csv"), stringsAsFactors=FALSE)
cis_dat <- read.csv(file.path(REV,"C_cis_TF_dat_500kb.csv"), stringsAsFactors=FALSE)
civw <- cis_mr[cis_mr$method=="Inverse variance weighted",]
add("cis_TF_denovo_2snps (Fig S9)", 2, paste(sort(cis_dat$SNP),collapse=";"),
    civw$b, civw$se, civw$pval, "IVW fixed-effect (de novo cis extraction, +/-500kb, r2<0.001)",
    "R1_C_cis_TF_MR.R")
# individual cis Wald ratios (Part 2 preview)
cis_dat$wald    <- cis_dat$beta.outcome / cis_dat$beta.exposure
cis_dat$wald_se <- cis_dat$se.outcome / abs(cis_dat$beta.exposure)
cis_dat$Fstat   <- (cis_dat$beta.exposure / cis_dat$se.exposure)^2
for (i in seq_len(nrow(cis_dat))) {
  add(sprintf("cis_%s_alone (Wald)", cis_dat$SNP[i]), 1, cis_dat$SNP[i],
      cis_dat$wald[i], cis_dat$wald_se[i],
      2*pnorm(-abs(cis_dat$wald[i]/cis_dat$wald_se[i])),
      sprintf("Wald ratio (F=%.0f)", cis_dat$Fstat[i]), "R1_C_cis_TF_MR.R")
}

# ---- MVMR rows (resolves 1.3) ----
mv <- read.csv(file.path(REV,"E_MVMR_extended.csv"), stringsAsFactors=FALSE)
add("MVMR_transferrin (2-way)", 10, "union transferrin+iron",
    mv$beta[mv$exposure=="transferrin"], mv$se[mv$exposure=="transferrin"],
    mv$pval[mv$exposure=="transferrin"], "MVMR-IVW t-based, df=8 (MVMR pkg)", "R1_E_MVMR_extended.R")
add("MVMR_serum_iron (2-way)", 10, "union transferrin+iron",
    mv$beta[mv$exposure=="serum_iron"], mv$se[mv$exposure=="serum_iron"],
    mv$pval[mv$exposure=="serum_iron"], "MVMR-IVW t-based, df=8 (MVMR pkg)", "R1_E_MVMR_extended.R")
eg <- read.csv(file.path(REV,"H_MVMR_Egger.csv"), stringsAsFactors=FALSE)
add("MVMR_Egger_transferrin_slope", 10, "union transferrin+iron",
    eg$estimate[eg$term=="transferrin"], eg$se[eg$term=="transferrin"],
    eg$pval[eg$term=="transferrin"], "MVMR-Egger", "R1_H_polish.R")

recon <- do.call(rbind, rows)
write.csv(recon, file.path(REV,"reconciliation.csv"), row.names=FALSE)

# ---- coloc posteriors (resolves 1.5) ----
co <- read.csv(file.path(REV,"D_coloc_TF_summary.csv"), stringsAsFactors=FALSE)
pp <- c(H0=co$PP.H0.abf, H1=co$PP.H1.abf, H2=co$PP.H2.abf, H3=co$PP.H3.abf, H4=co$PP.H4.abf)

# ---- rs8177240 Wald SE (resolves 1.4) ----
rs <- harm[harm$SNP=="rs8177240",]
rs_wald <- rs$beta.outcome/rs$beta.exposure
rs_wald_se <- rs$se.outcome/abs(rs$beta.exposure)

# ---- notes ----
FE <- function(nm) { r<-ivw(harm$wald[harm$SNP %in% subsets[[nm]]], harm$wald_se[harm$SNP %in% subsets[[nm]]]); r }
r6 <- FE("excl_HLA_FADS")

notes <- c(
"# Part 1 — reconciliation notes","",
sprintf("Generated %s. Source of truth for the primary panel: results/primary/harmonised_transferrin_tau_2026-03-28.csv (exposure ieu-a-1052 = Benyamin 2014 transferrin n=23,986; outcome ebi-a-GCST90095138 = Sarnowski 2022 total-tau n=14,721).", DATE_RUN),"",
"## 1.1 Cis-TF instrument count — RESOLVED: two different analyses, zero shared SNPs",
"- Fig S9 (Analysis C, script R1_C_cis_TF_MR.R): DE NOVO extraction of transferrin GWAS SNPs at p<5e-8 within +/-500 kb of the TF gene body (chr3:133,464,073-133,494,388), clumped r2<0.001 in EUR. -> 2 SNPs {rs3811658, rs17376530}, IVW beta=-0.060, SE=0.017, p=3.99e-4. Identical at +/-1 Mb.",
"- Fig S6 (Analysis B, script R1_B_leave_locus_out.R) row 'Cis-TF only (3 SNPs)': NOT a de novo extraction. It is the subset of the ORIGINAL 8-SNP primary panel hand-labelled TF/TFRC-region = {rs8177240, rs9990333, rs7646473}, IVW beta=-0.044, p=0.006.",
"- They share NO SNPs and are built differently. rs8177240 does not survive C's clumping because rs3811658 (849 bp away, high LD) tags the same signal and is retained instead.",
"- ACTION: report Analysis C as THE cis-only MR. Relabel the Fig S6 row 'TF/TFRC-region subset of primary panel (3 SNPs)', present as sensitivity, not a second cis definition.","",
"## 1.2 Excl. HLA + FADS (6 SNPs) — RESOLVED: fixed- vs random-effects IVW",
sprintf("- Point estimate identical: beta=%.4f across 6 SNPs {%s}.", r6$beta, paste(sort(subsets[['excl_HLA_FADS']]),collapse=', ')),
sprintf("- Fixed-effect IVW: SE=%.4f, p=%.4f (this is what Fig S6 / base-R Analysis B reports ~0.003).", r6$se_fe, r6$p_fe),
sprintf("- Multiplicative random-effect IVW: SE=%.4f, p=%.4f (this matches manuscript/Fig S1 ~0.006).", r6$se_re, r6$p_re),
sprintf("- Cause: residual heterogeneity (Q=%.2f, Q_p=%.3f, I2=%.0f%%) inflates the RE SE. Not a harmonization difference and not a mislabeled row.", r6$Q, r6$Q_p, r6$I2),
"- ACTION: standardize on ONE convention. The manuscript primary (all-8) uses TwoSampleMR default IVW = multiplicative random-effects; for all-8 there is no heterogeneity so FE=RE and both give 2.88e-3. For internal consistency, report the leave-locus-out table with RANDOM-EFFECTS IVW (the conservative, manuscript-consistent choice) -> excl_HLA_FADS p=0.006. Regenerate Fig S6 accordingly.","",
"## 1.3 Serum-iron MVMR p — RESOLVED: use 0.801 (t-based), retire 0.817 (z-based)",
sprintf("- 0.801 = MVMR package, t-distribution inference, df = n_snps - n_exposures = 10 - 2 = 8. SAME basis as the transferrin row (p=1.4e-3). Source: E_MVMR_extended.csv (serum_iron pval=%.4f).", mv$pval[mv$exposure=='serum_iron']),
"- 0.817 = MendelianRandomization::mr_mvivw z-based approximation (old Fig 3). 0.80 = 0.801 rounded to 2 dp.",
"- ACTION: report 0.801 (or 0.80) everywhere incl. Fig 3 panels A & B. It is not acceptable to show transferrin t-based (1.4e-3) alongside serum iron z-based (0.817) in the same figure.","",
"## 1.4 rs8177240 Wald SE — RESOLVED: 0.017",
sprintf("- From harmonized data: Wald beta = beta.outcome/beta.exposure = %.4f/%.4f = %.4f; SE = se.outcome/|beta.exposure| = %.4f/%.4f = %.5f -> 0.017 (3 dp). p=%.4f.",
        rs$beta.outcome, rs$beta.exposure, rs_wald, rs$se.outcome, abs(rs$beta.exposure), rs_wald_se, 2*pnorm(-abs(rs_wald/rs_wald_se))),
"- Fig S9 shows SE=0.016 (hardcoded 0.0165 that floors to 0.016). Manuscript's 0.017 is correct. ACTION: fix Fig S9 label to 0.017.","",
"## 1.5 Colocalization posteriors — RESOLVED: sum to 1.000",
sprintf("- PP.H0=%.3f  PP.H1=%.3f  PP.H2=%.3f  PP.H3=%.3f  PP.H4=%.3f  (sum=%.3f)", pp['H0'],pp['H1'],pp['H2'],pp['H3'],pp['H4'], sum(pp)),
"- Fig S5 currently omits H0 and H2. ACTION: show all five to 3 dp; note PP.H4=0.447 is below the conventional 0.8 threshold and the conditional H4/(H3+H4)=96.3% is reported because PP.H1=0.536 reflects limited tau-GWAS power, not evidence against colocalization.","",
"## Cross-cutting flag for Part 2",
sprintf("- Individual cis Wald ratios: rs3811658 beta=%.4f p=%.4f (F=%.0f); rs17376530 beta=%.4f p=%.4f (F=%.0f).",
        cis_dat$wald[cis_dat$SNP=='rs3811658'], 2*pnorm(-abs(cis_dat$wald[cis_dat$SNP=='rs3811658']/cis_dat$wald_se[cis_dat$SNP=='rs3811658'])), cis_dat$Fstat[cis_dat$SNP=='rs3811658'],
        cis_dat$wald[cis_dat$SNP=='rs17376530'], 2*pnorm(-abs(cis_dat$wald[cis_dat$SNP=='rs17376530']/cis_dat$wald_se[cis_dat$SNP=='rs17376530'])), cis_dat$Fstat[cis_dat$SNP=='rs17376530']),
"- The cis IVW is ~90% weighted on rs3811658 (which tags rs8177240); rs17376530 is directionally concordant but individually non-significant (p~0.10). LD r2 (Part 2, API) needed to quantify independence. Report plainly."
)
writeLines(notes, file.path(REV,"reconciliation_notes.md"))

cat("=== reconciliation.csv written:", nrow(recon), "rows ===\n\n")
print(recon[,c("analysis","n_snps","beta","se","pval","method")], row.names=FALSE)
cat("\n=== 1.2 excl_HLA_FADS: FE p=", signif(r6$p_fe,3), " RE p=", signif(r6$p_re,3),
    " (I2=", round(r6$I2), "%)\n", sep="")
cat("=== 1.5 coloc sum =", sum(pp), "\n")
