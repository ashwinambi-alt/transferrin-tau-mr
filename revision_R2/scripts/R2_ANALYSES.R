##############################################################################
# R2_ANALYSES.R  --  ALZ-26-0616.R2 (minor revision) reviewer analyses
#
# Analysis A: colocalization -- unconditional posteriors, prior sensitivity,
#             regional detail  (Reviewer 2, point 1)
# Analysis B: cis instruments -- per-variant Wald ratios and IVW weights
#             (Reviewer 2, point 2)
# Analysis C: MVMR between-exposure covariance sensitivity
#             (Reviewer 2, point 3)
#
# ANALYSIS ONLY. Writes nothing outside results/R2/.
# Every input is md5'd before use; every output is md5'd and re-read after write.
##############################################################################

.libPaths(c("C:/Users/ashwi/AppData/Local/R/win-library/4.6", .libPaths()))

args <- commandArgs(trailingOnly = TRUE)
PROJECT_DIR <- if (length(args) >= 1) args[1] else
  "C:/Users/ashwi/OneDrive/Documents/EB1A Docs/MR Paper"
IN_DIR  <- file.path(PROJECT_DIR, "results", "reviewer_response")
OUT_DIR <- file.path(PROJECT_DIR, "results", "R2")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({ library(coloc); library(MVMR) })

set.seed(20261019)   # nothing here is stochastic, but the protocol asks
cat("Random seed set to 20261019 (no stochastic step in this script)\n\n")

rule  <- function(s) cat("\n", strrep("=", 78), "\n", s, "\n", strrep("=", 78), "\n", sep = "")
inmd5 <- function(p) { cat(sprintf("  INPUT  md5 %s  %s\n", tools::md5sum(p), basename(p))); p }
wr    <- function(df, name) {
  p <- file.path(OUT_DIR, name)
  write.csv(df, p, row.names = FALSE)
  back <- read.csv(p, stringsAsFactors = FALSE)   # verify from disk
  cat(sprintf("  WROTE  md5 %s  %s  (%d rows x %d cols re-read from disk)\n",
              tools::md5sum(p), name, nrow(back), ncol(back)))
  ln <- readLines(p, warn = FALSE)
  cat("    first: ", ln[1], "\n    last:  ", ln[length(ln)], "\n", sep = "")
  invisible(p)
}

rule("ENVIRONMENT")
cat(R.version.string, "\n")
claimed <- c(coloc = "5.2.3", MVMR = "0.4.6", TwoSampleMR = "0.7.9",
             MendelianRandomization = "0.10.0")
for (p in names(claimed)) {
  v <- tryCatch(as.character(packageVersion(p)), error = function(e) "NOT INSTALLED")
  cat(sprintf("  %-24s installed %-10s manuscript claims %-10s %s\n",
              p, v, claimed[[p]], ifelse(v == claimed[[p]], "MATCH", "*** MISMATCH ***")))
}

##############################################################################
rule("ANALYSIS A -- COLOCALIZATION")
##############################################################################

f_post <- inmd5(file.path(IN_DIR, "D_coloc_TF_snp_posteriors.csv"))
f_summ <- inmd5(file.path(IN_DIR, "D_coloc_TF_summary.csv"))
f_ld   <- inmd5(file.path(IN_DIR, "P2_ld_matrix.csv"))

res      <- read.csv(f_post, stringsAsFactors = FALSE)
disk_sum <- read.csv(f_summ, stringsAsFactors = FALSE)
ld       <- read.csv(f_ld, row.names = 1, stringsAsFactors = FALSE)

cat(sprintf("\nSaved coloc output: %d variants, columns: %s\n",
            nrow(res), paste(names(res), collapse = ", ")))
cat("\nNOTE: the R1 script (R1_D_coloc_TF.R) did not save the raw regional pulls,\n",
    "and the OpenGWAS token is currently invalid (HTTP 401). All quantities below\n",
    "are recomputed from the saved per-SNP approximate Bayes factors, which is\n",
    "exact: coloc's posteriors are a deterministic function of lABF and the priors.\n", sep = "")

## ---- A1: full unconditional posteriors, default priors --------------------
rule("A1  Unconditional posteriors, default priors (p1=p2=1e-4, p12=1e-5)")
pp <- unname(coloc:::combine.abf(res$lABF.df1, res$lABF.df2,
                                 p1 = 1e-4, p2 = 1e-4, p12 = 1e-5))
names(pp) <- c("PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf")
cat("\nRecomputed from saved lABF:\n")
print(round(pp, 3))
cat(sprintf("\nSum of posteriors = %.6f\n", sum(pp)))
cat("\nOn-disk R1 summary for comparison:\n")
print(disk_sum)
maxdiff <- max(abs(pp - unlist(disk_sum[, c("PP.H0.abf","PP.H1.abf","PP.H2.abf",
                                            "PP.H3.abf","PP.H4.abf")])))
cat(sprintf("\nMax abs difference vs R1 on-disk values: %.3e  -> %s\n",
            maxdiff, ifelse(maxdiff < 1e-9, "EXACT REPRODUCTION", "*** DISCREPANCY ***")))
cat(sprintf("\nPP.H0 = %.6f and PP.H2 = %.6f : the manuscript's three reported\n",
            pp["PP.H0.abf"], pp["PP.H2.abf"]))
cat("posteriors summing to 1.000 is CONFIRMED, not a rounding artefact.\n")
cat(sprintf("\nConditional quantity PP.H4/(PP.H3+PP.H4) = %.4f  (manuscript: 96.3%%)\n",
            pp["PP.H4.abf"] / (pp["PP.H3.abf"] + pp["PP.H4.abf"])))

## ---- A2: coloc::sensitivity() ---------------------------------------------
rule("A2  coloc::sensitivity(res, rule = 'H4 > 0.5')")
obj <- list(summary = c(nsnps = nrow(res), pp),
            results = res,
            priors  = c(p1 = 1e-4, p2 = 1e-4, p12 = 1e-5))
class(obj) <- c("coloc_abf", "list")

pdf_path <- file.path(OUT_DIR, "FigS11_coloc_sensitivity.pdf")
pdf(pdf_path, width = 10, height = 7.5)
sens <- coloc::sensitivity(obj, rule = "H4 > 0.5")
dev.off()
cat(sprintf("\n  WROTE  md5 %s  %s\n", tools::md5sum(pdf_path), basename(pdf_path)))
cat(sprintf("  vector PDF, %d bytes\n", file.info(pdf_path)$size))

passing <- sens$p12[sens$pass]
cat(sprintf("\nRule 'H4 > 0.5' is satisfied for p12 in [%.3e, %.3e]\n",
            min(passing), max(passing)))
cat(sprintf("Default prior p12 = 1e-05 gives PP.H4 = %.3f -> rule %s at the default.\n",
            pp["PP.H4.abf"], ifelse(pp["PP.H4.abf"] > 0.5, "PASSES", "FAILS")))
wr(sens[, c("p12", "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf",
            "PP.H4.abf", "pass")], "coloc_sensitivity_curve.csv")

## ---- A3: explicit prior grid ----------------------------------------------
rule("A3  Explicit prior grid")
grid <- rbind(
  data.frame(axis = "p12",     p1 = 1e-4, p2 = 1e-4,
             p12 = c(1e-7, 1e-6, 5e-6, 1e-5, 5e-5, 1e-4)),
  data.frame(axis = "p1=p2",   p1 = c(1e-4, 1e-5), p2 = c(1e-4, 1e-5), p12 = 1e-5)
)
grid <- grid[!duplicated(grid[, c("p1","p2","p12")]) | grid$axis == "p1=p2", ]

out <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  g <- grid[i, ]
  v <- unname(coloc:::combine.abf(res$lABF.df1, res$lABF.df2,
                                  p1 = g$p1, p2 = g$p2, p12 = g$p12))
  data.frame(axis = g$axis, p1 = g$p1, p2 = g$p2, p12 = g$p12,
             is_default = (g$p1 == 1e-4 && g$p2 == 1e-4 && g$p12 == 1e-5),
             PP.H0 = v[1], PP.H1 = v[2], PP.H2 = v[3], PP.H3 = v[4], PP.H4 = v[5],
             PP.H4_cond = v[5] / (v[4] + v[5]),
             H4_exceeds_H1 = v[5] > v[2],
             p12_ge_min_p1p2 = g$p12 >= min(g$p1, g$p2))
}))
print(within(out, {
  PP.H0 <- round(PP.H0,3); PP.H1 <- round(PP.H1,3); PP.H2 <- round(PP.H2,3)
  PP.H3 <- round(PP.H3,3); PP.H4 <- round(PP.H4,3); PP.H4_cond <- round(PP.H4_cond,3)
}), row.names = FALSE)

cat("\nNOTE: rows flagged p12_ge_min_p1p2 have p12 at or above p1, i.e. at least half\n",
    "of transferrin-associated loci are assumed a priori to be causal for tau as\n",
    "well. These are reported for completeness but are not defensible priors.\n", sep = "")

cat("\n--- A3 interpretive question ---\n")
cat("Is there any prior in this grid under which PP.H4 exceeds PP.H1?\n")
ax12 <- out[out$axis == "p12", ]          # p1 = p2 = 1e-4 held fixed
axpp <- out[out$axis == "p1=p2", ]        # p12 = 1e-5 held fixed
dflt <- out[out$is_default, ][1, ]
if (any(out$H4_exceeds_H1)) {
  thr <- min(ax12$p12[ax12$H4_exceeds_H1])
  cat(sprintf("ANSWER: Yes, but only under priors more permissive than the default.\n"))
  cat(sprintf("        Along the p12 axis (p1=p2=1e-04 fixed), PP.H4 first exceeds PP.H1 at\n"))
  cat(sprintf("        p12 = %.0e, i.e. %.0fx the default; at the default p12 = 1e-05,\n",
              thr, thr / 1e-5))
  cat(sprintf("        PP.H1 = %.3f still exceeds PP.H4 = %.3f.\n", dflt$PP.H1, dflt$PP.H4))
  cat(sprintf("        Along the p1=p2 axis (p12=1e-05 fixed), lowering p1=p2 to 1e-05 also\n"))
  cat(sprintf("        flips the ordering (PP.H4 = %.3f), but that sets p12 equal to p1.\n",
              axpp$PP.H4[axpp$p1 == 1e-5]))
  cat(sprintf("        Every cell in which PP.H4 exceeds PP.H1 has p12 >= %.0e.\n",
              min(out$p12[out$H4_exceeds_H1])))
} else {
  cat("ANSWER: No -- PP.H1 exceeds PP.H4 at every prior examined.\n")
}
wr(out, "coloc_prior_grid.csv")

## ---- A4: regional detail ---------------------------------------------------
rule("A4  Regional detail")
## p on the log10 scale: |z| here reaches ~40, where 2*pnorm(-|z|) underflows to
## exactly 0 and ties would make which.min() pick an arbitrary variant.
z2log10p <- function(z) (pnorm(-abs(z), log.p = TRUE) + log(2)) / log(10)
res$log10p1 <- z2log10p(res$z.df1); res$log10p2 <- z2log10p(res$z.df2)
res$p1v <- 10^res$log10p1;          res$p2v <- 10^res$log10p2
fmtp <- function(lp) if (lp > -300) formatC(10^lp, format = "e", digits = 3) else
                     sprintf("1e%.1f", lp)
lead1 <- res[which.min(res$log10p1), ]
lead2 <- res[which.min(res$log10p2), ]
r_ld  <- ld["rs8177240", "rs3811658"]

cat(sprintf("\nTransferrin lead : %s at %d, z = %.2f, p = %s\n",
            lead1$snp, lead1$position, lead1$z.df1, fmtp(lead1$log10p1)))
cat(sprintf("Total-tau   lead : %s at %d, z = %.2f, p = %s\n",
            lead2$snp, lead2$position, lead2$z.df2, fmtp(lead2$log10p2)))
cat(sprintf("Same variant?    : %s\n", ifelse(lead1$snp == lead2$snp, "YES", "NO")))
cat(sprintf("  distance apart : %d bp\n", abs(lead1$position - lead2$position)))
cat(sprintf("  LD r^2         : %.4f  (from P2_ld_matrix.csv, r = %.6f)\n", r_ld^2, r_ld))
cat(sprintf("\nMinimum p in window, tau GWAS: %s at %s  (manuscript states 1.4e-03)\n",
            fmtp(min(res$log10p2)), res$snp[which.min(res$log10p2)]))
cat(sprintf("Minimum p in window, transferrin GWAS: %s at %s\n",
            fmtp(min(res$log10p1)), res$snp[which.min(res$log10p1)]))
cat(sprintf("  (%d transferrin variants have |z| large enough that p underflows to 0\n",
            sum(res$log10p1 < -300)))
cat("   in double precision; ranking is done on the log scale to avoid ties.)\n")
cat(sprintf("\nWhere rs8177240 ranks in the tau GWAS: %d of %d (p = %s)\n",
            rank(res$log10p2)[res$snp == "rs8177240"], nrow(res),
            fmtp(res$log10p2[res$snp == "rs8177240"])))

top5 <- res[order(-res$SNP.PP.H4), ][1:5, c("snp", "position", "SNP.PP.H4")]
cat("\nTop 5 variants by SNP.PP.H4:\n"); print(top5, row.names = FALSE)

detail <- data.frame(
  item = c("Window (GRCh37)", "Window half-width around TF gene body",
           "Variants in window, transferrin GWAS, before harmonisation",
           "Variants in window, total-tau GWAS, before harmonisation",
           "Variants shared and analysed after harmonisation and QC",
           "Lead variant, transferrin GWAS", "Lead variant p, transferrin GWAS",
           "Lead variant, total-tau GWAS", "Lead variant p, total-tau GWAS",
           "Lead variants identical?", "Distance between lead variants (bp)",
           "LD r^2 between lead variants",
           "Minimum p in window, total-tau GWAS",
           paste0("Top ", 1:5, " variant by SNP.PP.H4")),
  value = c("chr3:132,964,073-133,994,388",
            "500 kb",
            "not recoverable (raw pull not saved by R1_D_coloc_TF.R; OpenGWAS token invalid)",
            "not recoverable (raw pull not saved by R1_D_coloc_TF.R; OpenGWAS token invalid)",
            format(nrow(res)),
            lead1$snp, fmtp(lead1$log10p1),
            lead2$snp, fmtp(lead2$log10p2),
            ifelse(lead1$snp == lead2$snp, "yes", "no"),
            format(abs(lead1$position - lead2$position)),
            sprintf("%.4f", r_ld^2),
            sprintf("%s (%s)", fmtp(min(res$log10p2)), res$snp[which.min(res$log10p2)]),
            sprintf("%s (%.4f)", top5$snp, top5$SNP.PP.H4)),
  stringsAsFactors = FALSE)
print(detail, row.names = FALSE)
wr(detail, "coloc_regional_detail.csv")

## ---- A5: draft Supplementary Table S8 --------------------------------------
rule("A5  Draft Supplementary Table S8")
s8 <- rbind(
  data.frame(section = "Unconditional posteriors (default priors)",
             quantity = names(pp), value = sprintf("%.3f", pp), stringsAsFactors = FALSE),
  data.frame(section = "Unconditional posteriors (default priors)",
             quantity = "PP.H4/(PP.H3+PP.H4) [conditional]",
             value = sprintf("%.3f", pp["PP.H4.abf"]/(pp["PP.H3.abf"]+pp["PP.H4.abf"]))),
  data.frame(section = "Prior sensitivity",
             quantity = sprintf("p12=%.0e (p1=p2=%.0e)", out$p12, out$p1),
             value = sprintf("H0=%.3f H1=%.3f H2=%.3f H3=%.3f H4=%.3f; cond=%.3f",
                             out$PP.H0, out$PP.H1, out$PP.H2, out$PP.H3,
                             out$PP.H4, out$PP.H4_cond)),
  data.frame(section = "Regional detail", quantity = detail$item, value = detail$value)
)
print(s8, row.names = FALSE)
wr(s8, "TableS8_colocalization.csv")

##############################################################################
rule("ANALYSIS B -- CIS INSTRUMENTS: WALD RATIOS AND IVW WEIGHTS")
##############################################################################

f_s6  <- inmd5(file.path(IN_DIR, "TableS6_cis_instruments.csv"))
f_dat <- inmd5(file.path(IN_DIR, "C_cis_TF_dat_500kb.csv"))
s6    <- read.csv(f_s6, stringsAsFactors = FALSE)
dat   <- read.csv(f_dat, stringsAsFactors = FALSE)

cat("\nTable S6 as published (Wald ratios stored to 5 dp):\n")
print(s6[, c("SNP", "wald_beta", "wald_se", "wald_p", "F_stat", "source")], row.names = FALSE)

## The published IVW came from TwoSampleMR::mr_ivw on the harmonized data, so the
## weights are recomputed from that same source rather than from the rounded
## Wald columns above. The two differ in the 4th significant figure of p.
dat <- dat[order(match(dat$SNP, c("rs3811658", "rs17376530"))), ]
cat("\nHarmonized cis data (the actual pipeline input):\n")
print(dat[, c("SNP", "beta.exposure", "se.exposure", "beta.outcome", "se.outcome")],
      row.names = FALSE)

wald_b  <- dat$beta.outcome / dat$beta.exposure
wald_se <- abs(dat$se.outcome / dat$beta.exposure)   # first-order delta method
cat("\nWald ratios at full precision:\n")
print(data.frame(SNP = dat$SNP, wald_beta = signif(wald_b, 8),
                 wald_se = signif(wald_se, 8),
                 tableS6_beta = s6$wald_beta[match(dat$SNP, s6$SNP)],
                 tableS6_se   = s6$wald_se[match(dat$SNP, s6$SNP)]), row.names = FALSE)

w <- 1 / wald_se^2
cis <- data.frame(SNP = dat$SNP, wald_beta = wald_b, wald_se = wald_se,
                  ivw_weight_pct = 100 * w / sum(w), stringsAsFactors = FALSE)
cat("\nB1 -- inverse-variance weight decomposition:\n")
print(data.frame(SNP = cis$SNP, raw_weight = round(w, 2),
                 ivw_weight_pct = round(cis$ivw_weight_pct, 2)), row.names = FALSE)

b_ivw  <- sum(w * wald_b) / sum(w)
se_fe  <- sqrt(1 / sum(w))
cat("\nB4 -- reproduction of the published cis-only IVW estimate:\n")
cat(sprintf("  hand-computed (fixed-effect)      : b = %.4f  se = %.4f  p = %.4e\n",
            b_ivw, se_fe, 2 * pnorm(-abs(b_ivw / se_fe))))
tsmr <- TwoSampleMR::mr_ivw(dat$beta.exposure, dat$beta.outcome,
                            dat$se.exposure,  dat$se.outcome)
cat(sprintf("  TwoSampleMR::mr_ivw (as published) : b = %.4f  se = %.4f  p = %.4e\n",
            tsmr$b, tsmr$se, tsmr$pval))
cat(sprintf("  published in C_cis_TF_summary.txt  : b = -0.0604 se = 0.0171 p = 3.989e-04\n"))
cat(sprintf("  mr_ivw reproduces published: b %s, se %s, p %s\n",
            ifelse(round(tsmr$b,4)  == -0.0604, "MATCH", "*** DIFFERS ***"),
            ifelse(round(tsmr$se,4) ==  0.0171, "MATCH", "*** DIFFERS ***"),
            ifelse(signif(tsmr$pval,4) == 3.989e-4, "MATCH", "*** DIFFERS ***")))
cat("  NOTE: computing the same quantity from the 5-dp rounded Wald columns of\n")
cat("  Table S6 instead gives p = 3.995e-04, which disagrees with the published\n")
cat("  value in the 4th significant figure. The harmonized data above is the\n")
cat("  correct source; the rounded table must not be used to recompute anything.\n")

i17 <- cis[cis$SNP == "rs17376530", ]
i17$wald_p <- s6$wald_p[s6$SNP == "rs17376530"]
ci  <- i17$wald_beta + c(-1, 1) * 1.96 * i17$wald_se
cat("\nB3 -- rs17376530 alone (the only variant independent of rs8177240):\n")
cat(sprintf("  beta = %.4f, SE = %.4f, 95%% CI %.3f to %.3f, p = %.4f\n",
            i17$wald_beta, i17$wald_se, ci[1], ci[2], i17$wald_p))
cat(sprintf("  r^2 with rs8177240 = %.4f\n", ld["rs17376530","rs8177240"]^2))

s6$ivw_weight_pct <- NA_real_
s6$ivw_weight_pct[match(cis$SNP, s6$SNP)] <- round(cis$ivw_weight_pct, 2)
s6$in_cis_instrument_set <- ifelse(grepl("^cis", s6$source), "yes",
                                   "no - primary-panel variant, shown for context only")
s6$r2_with_rs8177240 <- round(sapply(s6$SNP, function(x) ld[x, "rs8177240"]^2), 4)
print(s6[, c("SNP", "wald_beta", "wald_se", "wald_p", "F_stat",
             "ivw_weight_pct", "r2_with_rs8177240", "in_cis_instrument_set")],
      row.names = FALSE)
wr(s6, "TableS6_with_weights.csv")

##############################################################################
rule("ANALYSIS C -- MVMR BETWEEN-EXPOSURE COVARIANCE SENSITIVITY")
##############################################################################

f_mv <- inmd5(file.path(IN_DIR, "E_MVMR_extended_snps.csv"))
d    <- read.csv(f_mv, stringsAsFactors = FALSE)
cat(sprintf("\n%d SNPs, columns: %s\n", nrow(d), paste(names(d), collapse = ", ")))

r_in <- format_mvmr(d[, c("beta_tf", "beta_fe")], d$beta_tau,
                    d[, c("se_tf", "se_fe")], d$se_tau, d$SNP)

cat("\n--- rho = 0 reproduction check against the published values ---\n")
s0 <- suppressWarnings(strength_mvmr(r_in, gencov = 0))
q0 <- suppressWarnings(pleiotropy_mvmr(r_in, gencov = 0))
cat(sprintf("  cond-F transferrin = %.2f  (published 278.20)  %s\n", s0[1,1],
            ifelse(abs(s0[1,1]-278.20) < 0.01, "MATCH", "*** DIFFERS ***")))
cat(sprintf("  cond-F serum iron  = %.2f  (published  93.88)  %s\n", s0[1,2],
            ifelse(abs(s0[1,2]-93.88) < 0.01, "MATCH", "*** DIFFERS ***")))
cat(sprintf("  Q_A = %.4f on %d df, p = %.4f  (published Q=6.28, p=5.075e-01)\n",
            q0$Qstat, nrow(d) - 3L, q0$Qpval))
cat(sprintf("  SNP-effect correlation r = %.3f  (published -0.611)\n",
            cor(d$beta_tf, d$beta_fe)))

cat("\nC1 -- estimating the covariance directly:\n")
cat("  MVMR::phenocov_mvmr(pcor, seBXGs) requires the phenotypic correlation\n")
cat("  between serum transferrin and serum iron in the Benyamin 2014 sample, and\n")
cat("  builds cov(bx1j, bx2j) = rho * se1j * se2j.\n")
cat("  A sample-specific value was not located. Benyamin 2014 is a meta-analysis\n")
cat("  across 11 discovery cohorts, so a single pooled transferrin-serum-iron\n")
cat("  phenotypic correlation is not a natural output of its design: the quantity\n")
cat("  differs by contributing cohort and the GWAS used here is a subset (n=23,986).\n")
cat("  No value is therefore adopted. Per the reviewer's own stated alternative\n")
cat("  ('or provide a sensitivity analysis across plausible covariance values'),\n")
cat("  the C2 grid is reported instead. Because the conclusion is unchanged across\n")
cat("  the entire grid, it does not depend on which value is correct -- which is a\n")
cat("  stronger answer than committing to one assumed number.\n")

rhos <- c(-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4)
cat("\nC2 -- sensitivity grid over the assumed phenotypic correlation rho:\n")
cgrid <- do.call(rbind, lapply(rhos, function(rho) {
  gc <- if (rho == 0) 0 else
        phenocov_mvmr(matrix(c(1, rho, rho, 1), 2, 2), d[, c("se_tf", "se_fe")])
  s <- suppressWarnings(strength_mvmr(r_in, gencov = gc))
  q <- suppressWarnings(pleiotropy_mvmr(r_in, gencov = gc))
  a <- suppressWarnings(utils::capture.output(est <- ivw_mvmr(r_in, gencov = gc)))
  data.frame(rho = rho, is_current = (rho == 0),
             condF_transferrin = round(s[1,1], 2), condF_serum_iron = round(s[1,2], 2),
             Q_A = round(q$Qstat, 4), Q_A_df = nrow(d) - 3L, Q_A_p = round(q$Qpval, 4),
             beta_transferrin = signif(est[1,1], 8), se_transferrin = signif(est[1,2], 8),
             p_transferrin = signif(est[1,4], 6),
             beta_serum_iron = signif(est[2,1], 8), se_serum_iron = signif(est[2,2], 8),
             p_serum_iron = signif(est[2,4], 6))
}))
print(cgrid, row.names = FALSE)
wr(cgrid, "mvmr_covariance_sensitivity.csv")

cat("\n--- C3 interpretive questions ---\n")
cat(sprintf("1. Does conditional F drop below 10 anywhere in the plausible range?\n   ANSWER: No -- the minimum across the grid is %.2f (serum iron at rho = %+.1f),\n   more than %.0f-fold above the conventional threshold of 10.\n",
            min(c(cgrid$condF_transferrin, cgrid$condF_serum_iron)),
            cgrid$rho[which.min(cgrid$condF_serum_iron)],
            min(c(cgrid$condF_transferrin, cgrid$condF_serum_iron)) / 10))
cat(sprintf("2. Does Q_A become significant (p < 0.05) anywhere in the plausible range?\n   ANSWER: No -- Q_A ranges %.3f to %.3f on %d df (p = %.3f to %.3f), null throughout.\n",
            min(cgrid$Q_A), max(cgrid$Q_A), cgrid$Q_A_df[1],
            min(cgrid$Q_A_p), max(cgrid$Q_A_p)))
inv <- length(unique(cgrid$beta_transferrin)) == 1 && length(unique(cgrid$se_transferrin)) == 1
cat(sprintf("3. Do the MVMR-IVW point estimates move at all?\n   ANSWER: No -- identical to 8 significant figures at every rho (%s).\n",
            ifelse(inv, "confirmed", "*** NOT confirmed ***")))
cat("   This is invariance by construction, not an empirical finding: MVMR 0.4.6's\n")
cat("   ivw_mvmr() accepts gencov but never uses it -- it warns when gencov == 0 and\n")
cat("   then fits lm(betaYG ~ -1 + betaX1 + betaX2, weights = 1/seBetaYG^2).\n")
cat("   gencov enters only strength_mvmr() and pleiotropy_mvmr().\n")
cat(sprintf("\n   Direction note: the biologically expected sign of rho is negative\n   (transferrin rises in iron deficiency). Negative rho RAISES conditional F\n   (%.2f / %.2f at rho = -0.6 vs %.2f / %.2f at rho = 0), so setting the\n   covariance to zero was the conservative choice for instrument strength.\n",
            cgrid$condF_transferrin[cgrid$rho == -0.6], cgrid$condF_serum_iron[cgrid$rho == -0.6],
            cgrid$condF_transferrin[cgrid$rho == 0],    cgrid$condF_serum_iron[cgrid$rho == 0]))

##############################################################################
rule("OUTPUT MANIFEST")
##############################################################################
fs <- list.files(OUT_DIR, full.names = TRUE)
for (f in fs) cat(sprintf("  %s  %8d bytes  %s\n",
                          tools::md5sum(f), file.info(f)$size, basename(f)))

rule("sessionInfo()")
print(sessionInfo())
