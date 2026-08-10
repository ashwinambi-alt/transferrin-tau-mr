##############################################################################
# R1_AUDIT.R — Full audit of reviewer-response analyses and manuscript numbers
#
# Four parts:
#   A. Reproduce manuscript headline numbers from cached RData
#   B. GWAS ID staleness check via ieugwasr::gwasinfo
#   C. Cross-check every number in SUMMARY_ALL_ANALYSES.md against source CSVs
#   D. Independent from-scratch re-computation of Analyses A and B
#
# Outputs:
#   results/reviewer_response/AUDIT_REPORT.md         (human-readable)
#   results/reviewer_response/AUDIT_details.csv       (structured findings)
##############################################################################

args <- commandArgs(trailingOnly = TRUE)
PROJECT_DIR <- if (length(args) >= 1) args[1] else normalizePath("..", mustWork = TRUE)
OUT_DIR     <- file.path(PROJECT_DIR, "results", "reviewer_response")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(ieugwasr); library(dplyr)
})

findings <- list()
add <- function(section, check, expected, observed, status, note = "") {
  if (length(observed) == 0) { observed <- NA; status <- "FAIL"; note <- paste("empty observed;", note) }
  if (length(expected) == 0) expected <- NA
  findings[[length(findings) + 1]] <<- data.frame(
    section = section, check = check,
    expected = as.character(expected[1]),
    observed = as.character(observed[1]),
    status = status, note = note,
    stringsAsFactors = FALSE
  )
}

TOL <- function(a, b, rel = 0.01, abs = 1e-4) {
  if (length(a) == 0 || length(b) == 0) return("FAIL")
  if (is.na(a) || is.na(b)) return("FAIL")
  if (is.logical(a) && is.logical(b)) return(if (identical(a, b)) "PASS" else "FAIL")
  if (is.character(a) || is.character(b)) return(if (identical(as.character(a), as.character(b))) "PASS" else "FAIL")
  a <- as.numeric(a); b <- as.numeric(b)
  ok <- (base::abs(a - b) < abs) ||
        (base::abs(a - b) / max(base::abs(a), base::abs(b), 1e-12) < rel)
  if (isTRUE(ok)) "PASS" else "FAIL"
}

# ============================================================
# A. Reproduce manuscript headline numbers from cached RData
# ============================================================
cat("== A. Reproduce manuscript headlines from cached RData ==\n")
rdata <- file.path(PROJECT_DIR, "data", "FINAL_analysis_results.RData")
if (file.exists(rdata)) {
  e <- new.env(); load(rdata, envir = e)
  cat("  loaded objects:", paste(ls(e), collapse = ", "), "\n")

  # A1. Primary IVW transferrin -> total tau
  pr <- e$primary_results
  ivw <- pr[pr$method == "Inverse variance weighted" &
             pr$pair == "transferrin_circ_total_tau", ]
  add("A", "Primary IVW β (manuscript = −0.039)",
      -0.039, signif(ivw$b, 4), TOL(ivw$b, -0.039, abs = 0.001))
  add("A", "Primary IVW SE (manuscript = 0.013)",
      0.013, signif(ivw$se, 4), TOL(ivw$se, 0.013, abs = 0.001))
  add("A", "Primary IVW p (manuscript = 2.88e-3)",
      2.88e-3, signif(ivw$pval, 4), TOL(ivw$pval, 2.88e-3, rel = 0.05))
  add("A", "Primary IVW nsnp (manuscript = 8)",
      8, ivw$nsnp, TOL(ivw$nsnp, 8, abs = 0))

  # A2. MVMR from cached
  mvmr_csv <- file.path(PROJECT_DIR, "results", "mvmr", "FINAL_mvmr_2026-03-28.csv")
  if (file.exists(mvmr_csv)) {
    m <- read.csv(mvmr_csv)
    tf_mv <- m[m$exposure == "Transferrin", ]
    fe_mv <- m[m$exposure == "Serum Iron", ]
    add("A", "MVMR transferrin β (manuscript = −0.044)",
        -0.044, signif(tf_mv$beta, 4), TOL(tf_mv$beta, -0.044, abs = 0.002))
    add("A", "MVMR transferrin p (manuscript = 2.29e-5)",
        2.29e-5, signif(tf_mv$pval, 4), TOL(tf_mv$pval, 2.29e-5, rel = 0.05))
    add("A", "MVMR serum iron p (manuscript = 0.817)",
        0.817, signif(fe_mv$pval, 4), TOL(fe_mv$pval, 0.817, abs = 0.01))
  } else {
    add("A", "MVMR cached CSV present", "exists", "missing", "FAIL")
  }

  # A3. Sensitivity results
  sens <- e$sensitivity_results
  add("A", "Egger intercept p (manuscript = 0.202)",
      0.202, signif(sens$egger_intercept_p, 3), TOL(sens$egger_intercept_p, 0.202, abs = 0.02))
  add("A", "Cochran Q p (manuscript = 0.544)",
      0.544, signif(sens$cochran_q_p, 3), TOL(sens$cochran_q_p, 0.544, abs = 0.02))
  add("A", "I² (manuscript = 0%)",
      0, signif(sens$I2, 3), TOL(sens$I2, 0, abs = 0.5))
  add("A", "MR-PRESSO global p (manuscript = 0.57)",
      0.57, signif(sens$presso_global_p, 3), TOL(sens$presso_global_p, 0.57, abs = 0.05))
  add("A", "Steiger correct SNPs (manuscript = 8/8)",
      "8/8",
      sprintf("%d/%d", sens$steiger_correct, sens$steiger_total),
      TOL(sens$steiger_correct, 8, abs = 0))
  add("A", "LOO removals losing significance (manuscript = 1 (rs8177240))",
      1, sens$loo_lose_sig, TOL(sens$loo_lose_sig, 1, abs = 0))

  # A4. Weighted median (from primary_results)
  wm <- pr[pr$method == "Weighted median" & pr$pair == "transferrin_circ_total_tau", ]
  add("A", "Weighted median p (manuscript = 0.001)",
      1.43e-3, signif(wm$pval, 4), TOL(wm$pval, 1.43e-3, rel = 0.1))

  # A5. Controls
  ctrl <- e$controls
  if (!is.null(ctrl$positive)) {
    add("A", "Positive control iron→Hb β (manuscript = +0.704)",
        0.704, signif(ctrl$positive$b, 4),
        TOL(ctrl$positive$b, 0.704, abs = 0.02))
  }
  if (!is.null(ctrl$negative)) {
    add("A", "Negative control iron→height p (manuscript = 0.191)",
        0.191, signif(ctrl$negative$pval, 3),
        TOL(ctrl$negative$pval, 0.191, abs = 0.03))
  }
} else {
  add("A", "cached RData exists", "exists", "missing", "FAIL",
      "cannot reproduce manuscript numbers")
}

# ============================================================
# B. GWAS ID staleness check
# ============================================================
cat("\n== B. GWAS ID staleness check ==\n")
gwas_ids <- c(
  "ieu-a-1049"         = "Benyamin 2014 serum iron",
  "ieu-a-1050"         = "Benyamin 2014 ferritin",
  "ieu-a-1051"         = "Benyamin 2014 TSAT",
  "ieu-a-1052"         = "Benyamin 2014 transferrin",
  "ieu-b-5115"         = "Bell 2021 ferritin (large)",
  "ebi-a-GCST90095138" = "Sarnowski 2022 total tau",
  "ebi-a-GCST90028995" = "Vuckovic 2020 haemoglobin (pos ctrl)",
  "ebi-a-GCST90029008" = "Vuckovic 2020 (neg ctrl)"
)

# Manuscript-cited sample sizes
expected_n <- c(
  "ieu-a-1049" = 23986, "ieu-a-1050" = 23986,
  "ieu-a-1051" = 23986, "ieu-a-1052" = 23986,
  "ieu-b-5115" = 246139,
  "ebi-a-GCST90095138" = 14721
)

gwas_info <- tryCatch(gwasinfo(names(gwas_ids)),
                       error = function(e) { cat("  gwasinfo err:",
                         conditionMessage(e), "\n"); NULL })
if (!is.null(gwas_info)) {
  print(gwas_info[, intersect(c("id","trait","year","author","sample_size",
                                 "population","nsnp"), names(gwas_info))])
  for (id in names(gwas_ids)) {
    row <- gwas_info[gwas_info$id == id, ]
    if (nrow(row) == 0) {
      add("B", paste0("GWAS ", id, " (", gwas_ids[id], ") resolvable"),
          "resolves", "MISSING", "FAIL",
          "OpenGWAS no longer serves this ID — dataset may have been withdrawn")
      next
    }
    add("B", paste0(id, " trait"),
        gwas_ids[id], row$trait, "INFO")
    if (id %in% names(expected_n)) {
      add("B", paste0(id, " sample size (manuscript = ", expected_n[id], ")"),
          expected_n[id], row$sample_size,
          if (abs(as.numeric(row$sample_size) - expected_n[id]) /
              expected_n[id] < 0.02) "PASS" else "WARN",
          "OpenGWAS metadata may round; ±2% considered OK")
    }
  }
} else {
  add("B", "gwasinfo call succeeded", "yes", "no", "FAIL")
}

# ============================================================
# C. Cross-check summary values against source CSVs
# ============================================================
cat("\n== C. Cross-check reviewer-response summary numbers ==\n")

check_csv <- function(path, checks) {
  if (!file.exists(path)) {
    add("C", paste("file:", basename(path)), "exists", "missing", "FAIL")
    return(NULL)
  }
  d <- read.csv(path, stringsAsFactors = FALSE)
  for (chk in checks) chk(d, basename(path))
}

# C.A rs8177240 contribution CSV
check_csv(file.path(OUT_DIR, "A_rs8177240_contribution.csv"), list(function(d, f) {
  rr <- d[d$SNP == "rs8177240", ]
  add("C-A", paste0(f, ": rs8177240 weight_pct (report = 63.5)"),
      63.5, rr$weight_pct, TOL(rr$weight_pct, 63.5, abs = 0.1))
  rr2 <- d[d$SNP == "rs1800562", ]
  add("C-A", paste0(f, ": rs1800562 weight_pct (report = 27.15)"),
      27.15, rr2$weight_pct, TOL(rr2$weight_pct, 27.15, abs = 0.1))
  add("C-A", paste0(f, ": sum of weight_pct == 100"),
      100.0, sum(d$weight_pct), TOL(sum(d$weight_pct), 100, abs = 0.1))
}))

# C.B leave-locus-out
check_csv(file.path(OUT_DIR, "B_leave_locus_out.csv"), list(function(d, f) {
  joint <- d[grepl("joint", d$subset), ]
  add("C-B", paste0(f, ": joint HLA+FADS+HFE+NAT2 ivw_beta (report −0.043)"),
      -0.043, joint$ivw_beta, TOL(joint$ivw_beta, -0.043, abs = 0.002))
  add("C-B", paste0(f, ": joint ivw_pval (report 0.0073)"),
      0.0073, joint$ivw_pval, TOL(joint$ivw_pval, 0.0073, abs = 0.001))
  cis3 <- d[grepl("cis_TF_3snps", d$subset), ]
  add("C-B", paste0(f, ": cis_TF_3snps ivw_beta (report −0.044)"),
      -0.044, cis3$ivw_beta, TOL(cis3$ivw_beta, -0.044, abs = 0.002))
  cis3wm <- cis3$wm_pval
  add("C-B", paste0(f, ": cis_TF_3snps wm_pval (report 0.0037)"),
      0.0037, cis3wm, TOL(cis3wm, 0.0037, rel = 0.05))
}))

# C.C cis-TF results
check_csv(file.path(OUT_DIR, "C_cis_TF_mr_500kb.csv"), list(function(d, f) {
  ivw <- d[d$method == "Inverse variance weighted", ]
  add("C-C", paste0(f, ": cis-TF IVW β (report −0.060)"),
      -0.060, signif(ivw$b, 4), TOL(ivw$b, -0.060, abs = 0.002))
  add("C-C", paste0(f, ": cis-TF IVW p (report 3.99e-4)"),
      3.99e-4, signif(ivw$pval, 4), TOL(ivw$pval, 3.99e-4, rel = 0.05))
  add("C-C", paste0(f, ": cis-TF nsnp (report 2)"),
      2, ivw$nsnp, TOL(ivw$nsnp, 2, abs = 0))
}))

# C.D coloc
check_csv(file.path(OUT_DIR, "D_coloc_TF_summary.csv"), list(function(d, f) {
  add("C-D", paste0(f, ": PP.H4 (report 0.447)"),
      0.447, signif(d$PP.H4.abf, 3), TOL(d$PP.H4.abf, 0.447, abs = 0.005))
  add("C-D", paste0(f, ": PP.H3 (report 0.017)"),
      0.017, signif(d$PP.H3.abf, 3), TOL(d$PP.H3.abf, 0.017, abs = 0.005))
  add("C-D", paste0(f, ": PP.H1 (report 0.536)"),
      0.536, signif(d$PP.H1.abf, 3), TOL(d$PP.H1.abf, 0.536, abs = 0.005))
  add("C-D", paste0(f, ": sum PP ≈ 1"),
      1.0, sum(d[, grepl("^PP\\.", names(d))]),
      TOL(sum(d[, grepl("^PP\\.", names(d))]), 1.0, abs = 0.01))
}))

# C.E MVMR extended
check_csv(file.path(OUT_DIR, "E_MVMR_extended.csv"), list(function(d, f) {
  tf <- d[d$exposure == "transferrin", ]
  add("C-E", paste0(f, ": MVMR transferrin β (report −0.044)"),
      -0.044, signif(tf$beta, 4), TOL(tf$beta, -0.044, abs = 0.002))
  add("C-E", paste0(f, ": transferrin cond-F (report 278.2)"),
      278.2, signif(tf$cond_F, 4), TOL(tf$cond_F, 278.2, rel = 0.01))
  fe <- d[d$exposure == "serum_iron", ]
  add("C-E", paste0(f, ": serum iron cond-F (report 93.9)"),
      93.9, signif(fe$cond_F, 4), TOL(fe$cond_F, 93.9, rel = 0.01))
  add("C-E", paste0(f, ": serum iron p (report 0.80)"),
      0.80, signif(fe$pval, 4), TOL(fe$pval, 0.80, abs = 0.02))
}))

# C.F MVMR 4-way
check_csv(file.path(OUT_DIR, "F_MVMR_4way.csv"), list(function(d, f) {
  # MVMR pkg writes positional labels "exposure1", "exposure2", etc.
  # 3-way order was transferrin, ferritin, tsat (per Analysis F script)
  m3 <- d[d$model == "3-way (transferrin + ferritin + tsat)", ]
  m3tf <- m3[m3$exposure == "exposure1", ]  # transferrin is position 1
  add("C-F", paste0(f, ": 3-way transferrin β (report −0.046)"),
      -0.046, signif(m3tf$beta, 4), TOL(m3tf$beta, -0.046, abs = 0.002))
  add("C-F", paste0(f, ": 3-way transferrin p (report 0.014)"),
      0.014, signif(m3tf$pval, 4), TOL(m3tf$pval, 0.014, abs = 0.001))
  add("C-F", paste0(f, ": 3-way transferrin cond-F (report 22.7)"),
      22.7, signif(m3tf$cond_F, 4), TOL(m3tf$cond_F, 22.7, abs = 0.5))
  # ferritin (exposure2) and tsat (exposure3) should both be null
  m3fe <- m3[m3$exposure == "exposure2", ]  # ferritin
  m3ts <- m3[m3$exposure == "exposure3", ]  # tsat
  add("C-F", paste0(f, ": 3-way ferritin p (report ~0.82)"),
      0.824, signif(m3fe$pval, 3), TOL(m3fe$pval, 0.824, abs = 0.01))
  add("C-F", paste0(f, ": 3-way TSAT p (report ~0.83)"),
      0.833, signif(m3ts$pval, 3), TOL(m3ts$pval, 0.833, abs = 0.01))
}))

# C.G GWAS Catalog summary
check_csv(file.path(OUT_DIR, "G_gwas_catalog_summary.csv"), list(function(d, f) {
  tf_row <- d[d$snp == "rs8177240", ]
  add("C-G", paste0(f, ": rs8177240 n_iron (report 5)"),
      5, tf_row$n_iron, TOL(tf_row$n_iron, 5, abs = 0))
  hfe_row <- d[d$snp == "rs1800562", ]
  add("C-G", paste0(f, ": rs1800562 n_hits (report 183)"),
      183, hfe_row$n_hits, TOL(hfe_row$n_hits, 183, abs = 0))
  add("C-G", paste0(f, ": rs8177240 concern (report LOW)"),
      "LOW", tf_row$pleiotropy_concern,
      if (grepl("LOW", tf_row$pleiotropy_concern)) "PASS" else "FAIL")
}))

# C.H polish
check_csv(file.path(OUT_DIR, "H_MVMR_Egger.csv"), list(function(d, f) {
  int <- d[d$term == "intercept", ]
  add("C-H", paste0(f, ": MVMR-Egger intercept p (report 0.244)"),
      0.244, signif(int$pval, 3), TOL(int$pval, 0.244, abs = 0.005))
  tfsl <- d[d$term == "transferrin", ]
  add("C-H", paste0(f, ": MVMR-Egger transferrin slope β (report −0.056)"),
      -0.056, signif(tfsl$estimate, 4), TOL(tfsl$estimate, -0.056, abs = 0.002))
  add("C-H", paste0(f, ": transferrin slope p (report 1.65e-4)"),
      1.65e-4, signif(tfsl$pval, 4), TOL(tfsl$pval, 1.65e-4, rel = 0.05))
}))

check_csv(file.path(OUT_DIR, "H_interaction_tests.csv"), list(function(d, f) {
  apo <- d[grepl("APOE4 non-carriers vs carriers", d$contrast) & grepl("p-tau", d$contrast), ]
  add("C-H", paste0(f, ": APOE4 p-tau interaction p_Q (report 0.095)"),
      0.0952, signif(apo$p_Q, 3), TOL(apo$p_Q, 0.0952, abs = 0.005))
  ab <- d[grepl("APOE4 non-carriers vs carriers", d$contrast) & grepl("Abeta", d$contrast), ]
  add("C-H", paste0(f, ": Aβ42 APOE4 interaction p_Q (report 0.96)"),
      0.96, signif(ab$p_Q, 3), TOL(ab$p_Q, 0.96, abs = 0.005))
}))

# ============================================================
# D. Independent from-scratch re-computation of Analyses A and B
# ============================================================
cat("\n== D. From-scratch re-computation from harmonised CSV ==\n")
h <- read.csv(file.path(PROJECT_DIR, "results", "primary",
                        "harmonised_transferrin_tau_2026-03-28.csv"),
              stringsAsFactors = FALSE)
h <- h[h$mr_keep == TRUE, ]
add("D", "harmonised CSV row count",
    8, nrow(h), TOL(nrow(h), 8, abs = 0))

# D1. Independent IVW recompute
h$wald   <- h$beta.outcome / h$beta.exposure
h$wald_se <- h$se.outcome / abs(h$beta.exposure)
w <- 1 / h$wald_se^2
ivw_b <- sum(w * h$wald) / sum(w)
ivw_se <- sqrt(1 / sum(w))
ivw_p  <- 2 * pnorm(-abs(ivw_b / ivw_se))
add("D", "IVW β from independent recompute (should = −0.039)",
    -0.039, signif(ivw_b, 4), TOL(ivw_b, -0.039, abs = 0.001))
add("D", "IVW p from independent recompute (should = 2.88e-3)",
    2.88e-3, signif(ivw_p, 4), TOL(ivw_p, 2.88e-3, rel = 0.05))

# D2. Weight percent recompute
w_pct <- 100 * w / sum(w)
rs8_pct <- w_pct[h$SNP == "rs8177240"]
add("D", "rs8177240 IVW weight % recompute (should = 63.5)",
    63.5, signif(rs8_pct, 4), TOL(rs8_pct, 63.5, abs = 0.1))

# D3. LOO of rs8177240 recompute
sub <- h[h$SNP != "rs8177240", ]
ws <- 1 / (sub$se.outcome / abs(sub$beta.exposure))^2
subb <- sum(ws * (sub$beta.outcome / sub$beta.exposure)) / sum(ws)
subse <- sqrt(1 / sum(ws))
subp  <- 2 * pnorm(-abs(subb / subse))
add("D", "IVW β excl rs8177240 recompute (should = −0.020)",
    -0.020, signif(subb, 4), TOL(subb, -0.020, abs = 0.002))
add("D", "IVW p excl rs8177240 recompute (should = 0.358)",
    0.358, signif(subp, 3), TOL(subp, 0.358, abs = 0.01))

# D4. Allele-alignment audit — for each SNP, effect_allele.exposure should equal
# effect_allele.outcome after harmonise_data(action=2). This is the actual
# invariant (not the beta sign, which reflects the SNP's biology).
allele_consistent <- with(h,
  all(effect_allele.exposure == effect_allele.outcome))
add("D", "harmonise action=2: effect_allele.exposure == effect_allele.outcome per SNP",
    TRUE, allele_consistent, TOL(allele_consistent, TRUE, abs = 0))
# And under action=2, mixed +/− exposure betas are EXPECTED (they reflect
# which allele each SNP was reported against in the original GWAS)
n_neg <- sum(h$beta.exposure < 0)
add("D", "beta.exposure sign distribution",
    "3 negative (rs1800562, rs7646473, rs9990333)",
    sprintf("%d negative (%s)", n_neg,
            paste(h$SNP[h$beta.exposure < 0], collapse = ", ")),
    if (n_neg == 3) "PASS" else "WARN",
    "TwoSampleMR harmonise action=2 preserves the exposure GWAS's reported effect direction; mixed signs are normal.")

# D5. Cross-check that the sensitivity CSV "excl rs8177240" matches manuscript
sens_csv <- file.path(PROJECT_DIR, "results", "sensitivity",
                       "sensitivity_no_rs8177240_2026-03-28.csv")
if (file.exists(sens_csv)) {
  s <- read.csv(sens_csv)
  ivw_row <- s[grepl("Inverse variance weighted", s$method), ]
  add("D", "cached excl-rs8177240 β matches manuscript −0.020",
      -0.020, signif(ivw_row$b[1], 4),
      TOL(ivw_row$b[1], -0.020, abs = 0.002))
  add("D", "cached excl-rs8177240 p matches manuscript 0.358",
      0.358, signif(ivw_row$pval[1], 3),
      TOL(ivw_row$pval[1], 0.358, abs = 0.005))
}

# ============================================================
# Save findings
# ============================================================
findings_df <- do.call(rbind, findings)
write.csv(findings_df, file.path(OUT_DIR, "AUDIT_details.csv"), row.names = FALSE)

n_pass <- sum(findings_df$status == "PASS")
n_fail <- sum(findings_df$status == "FAIL")
n_warn <- sum(findings_df$status == "WARN")
n_info <- sum(findings_df$status == "INFO")

cat(sprintf("\nAudit results: %d PASS, %d FAIL, %d WARN, %d INFO\n",
            n_pass, n_fail, n_warn, n_info))
if (n_fail > 0) {
  cat("\n== FAILURES ==\n")
  print(findings_df[findings_df$status == "FAIL", ])
}
if (n_warn > 0) {
  cat("\n== WARNINGS ==\n")
  print(findings_df[findings_df$status == "WARN", ])
}

# Human-readable report
lines <- c(
  "# Audit Report — reviewer-response analyses",
  sprintf("Run: %s", Sys.time()),
  "",
  sprintf("Summary: %d PASS, %d FAIL, %d WARN, %d INFO", n_pass, n_fail, n_warn, n_info),
  "",
  "## Section A — Reproduce manuscript headline numbers",
  ""
)
tbl_sect <- function(sec) {
  s <- findings_df[grepl(paste0("^", sec), findings_df$section), ]
  if (nrow(s) == 0) return("(no checks in this section)")
  out <- c("| Check | Expected | Observed | Status |",
           "|---|---|---|---|")
  for (i in seq_len(nrow(s))) {
    out <- c(out, sprintf("| %s | %s | %s | %s |",
                          s$check[i], s$expected[i], s$observed[i], s$status[i]))
  }
  out
}
lines <- c(lines, tbl_sect("A"), "",
           "## Section B — GWAS ID staleness", "", tbl_sect("B"), "",
           "## Section C — Cross-check summary values against source CSVs", "",
           tbl_sect("C"), "",
           "## Section D — Independent from-scratch recompute", "",
           tbl_sect("D"))
writeLines(lines, file.path(OUT_DIR, "AUDIT_REPORT.md"))
cat(sprintf("\nWritten: %s\nWritten: %s\n",
            file.path(OUT_DIR, "AUDIT_REPORT.md"),
            file.path(OUT_DIR, "AUDIT_details.csv")))
