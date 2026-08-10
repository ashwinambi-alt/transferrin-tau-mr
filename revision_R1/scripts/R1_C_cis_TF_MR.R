##############################################################################
# R1_C_cis_TF_MR.R
#
# Reviewer response C: Cis-only Mendelian randomization at the TF locus.
#
# Extracts all transferrin-associated SNPs within +/- 500 kb of the TF gene
# body (chr3:133,464,073 - 133,494,388 GRCh37) at p < 5e-8, clumps them at
# r2 < 0.001, and runs IVW / weighted median / MR-Egger against circulating
# total-tau (ebi-a-GCST90095138).
#
# Cis-only MR is the cleanest test of transferrin-specific biology because
# any horizontal-pleiotropy signal from FADS / HLA / HFE / NAT2 is excluded
# by construction.
#
# Sensitivity: also runs ±1 Mb.
##############################################################################

args <- commandArgs(trailingOnly = TRUE)
PROJECT_DIR <- if (length(args) >= 1) args[1] else normalizePath("..", mustWork = TRUE)
OUT_DIR     <- file.path(PROJECT_DIR, "results", "reviewer_response")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(ieugwasr)
  library(TwoSampleMR)
  library(dplyr)
})

stopifnot(nchar(Sys.getenv("OPENGWAS_JWT")) > 100)

# TF gene body coordinates on GRCh37 (Ensembl)
TF_CHR   <- 3L
TF_START <- 133464073L
TF_END   <- 133494388L

run_window <- function(window_kb, tag) {
  win_start <- TF_START - window_kb * 1000
  win_end   <- TF_END   + window_kb * 1000
  cat(sprintf("\n== Cis-only TF MR (+/- %d kb, %s Mb window) ==\n",
              window_kb, tag))
  cat(sprintf("  Region: chr%d:%d-%d (%.2f Mb)\n",
              TF_CHR, win_start, win_end, (win_end - win_start) / 1e6))

  # ---- pull ALL genome-wide-significant transferrin SNPs, then filter to cis ----
  # extract_instruments already clumps; but we want to bypass the clumping so we
  # can clump only within cis. Use ieugwasr::tophits with clump=FALSE.
  cat("  Fetching transferrin SNPs at p<5e-8 without clumping...\n")
  tf_all <- tryCatch(
    tophits(id = "ieu-a-1052", pval = 5e-8, r2 = 1, kb = 0, clump = 0),
    error = function(e) { cat("    ERR:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(tf_all)) stop("Failed to fetch transferrin instruments")
  cat(sprintf("    Retrieved %d SNPs genome-wide\n", nrow(tf_all)))

  # Standardise column names (tophits uses OpenGWAS schema)
  # chr / position columns are usually 'chr' and 'position' in ieugwasr output
  chr_col <- intersect(c("chr", "chromosome"), names(tf_all))[1]
  pos_col <- intersect(c("position", "pos", "bp"), names(tf_all))[1]
  if (is.na(chr_col) || is.na(pos_col)) {
    print(head(tf_all)); stop("Cannot identify chr/pos columns")
  }
  tf_all[[chr_col]] <- as.integer(tf_all[[chr_col]])
  tf_all[[pos_col]] <- as.integer(tf_all[[pos_col]])

  cis <- tf_all[tf_all[[chr_col]] == TF_CHR &
                tf_all[[pos_col]] >= win_start &
                tf_all[[pos_col]] <= win_end, ]
  cat(sprintf("    %d SNPs fall in the cis window\n", nrow(cis)))
  if (nrow(cis) == 0) return(NULL)

  # ---- clump inside cis window ----
  cat("  Clumping at r2<0.001...\n")
  # ld_clump uses local reference by default; pass through OpenGWAS
  cis_clump <- tryCatch(
    ld_clump(dplyr::tibble(rsid = cis$rsid, pval = cis$p, id = cis$id),
             clump_r2 = 0.001, clump_kb = 10000,
             pop = "EUR", opengwas_jwt = Sys.getenv("OPENGWAS_JWT")),
    error = function(e) { cat("    LD clump error:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(cis_clump)) {
    cat("    Clumping failed — using unclumped set for reporting\n")
    cis_clump <- data.frame(rsid = cis$rsid)
  }
  cat(sprintf("    %d independent cis SNPs after clumping\n", nrow(cis_clump)))
  if (nrow(cis_clump) == 0) return(NULL)

  # ---- build exposure dataframe from the tophits pull ----
  cis_kept <- cis[cis$rsid %in% cis_clump$rsid, ]
  # Convert to TwoSampleMR exposure format
  # tophits columns typically: rsid, ea, nea, eaf, beta, se, p, n
  exp_df <- format_data(
    data.frame(
      SNP            = cis_kept$rsid,
      beta           = cis_kept$beta,
      se             = cis_kept$se,
      effect_allele  = cis_kept$ea,
      other_allele   = cis_kept$nea,
      eaf            = cis_kept$eaf,
      pval           = cis_kept$p,
      samplesize     = cis_kept$n,
      chr            = cis_kept[[chr_col]],
      pos            = cis_kept[[pos_col]],
      exposure       = "transferrin (cis-TF)"
    ), type = "exposure",
    snp_col="SNP", beta_col="beta", se_col="se",
    effect_allele_col="effect_allele", other_allele_col="other_allele",
    eaf_col="eaf", pval_col="pval", samplesize_col="samplesize",
    chr_col="chr", pos_col="pos"
  )
  exp_df$F_stat <- (exp_df$beta.exposure / exp_df$se.exposure)^2
  cat(sprintf("    F-stat range: %.1f – %.1f (mean %.1f)\n",
              min(exp_df$F_stat), max(exp_df$F_stat), mean(exp_df$F_stat)))

  # ---- pull outcome ----
  cat("  Extracting outcome (total-tau)...\n")
  out_df <- extract_outcome_data(snps = exp_df$SNP,
                                 outcomes = "ebi-a-GCST90095138")
  if (is.null(out_df) || nrow(out_df) == 0) {
    cat("    No outcome data — abort\n"); return(NULL)
  }
  cat(sprintf("    Retrieved %d outcome SNPs\n", nrow(out_df)))

  # ---- harmonise + MR ----
  dat <- harmonise_data(exp_df, out_df, action = 2)
  dat <- dat[dat$mr_keep == TRUE, ]
  cat(sprintf("    %d SNPs pass harmonisation\n", nrow(dat)))
  if (nrow(dat) < 1) return(NULL)

  # Multiple methods only meaningful with >=3 SNPs
  methods <- if (nrow(dat) >= 3) {
    c("mr_ivw", "mr_egger_regression", "mr_weighted_median",
      "mr_weighted_mode", "mr_simple_mode")
  } else if (nrow(dat) == 2) {
    c("mr_ivw")
  } else {
    c("mr_wald_ratio")
  }
  mr_res <- mr(dat, method_list = methods)
  print(mr_res)

  # Egger intercept + Q + PRESSO if >=3
  extras <- list()
  if (nrow(dat) >= 3) {
    extras$egger_intercept <- mr_pleiotropy_test(dat)
    extras$Q               <- mr_heterogeneity(dat)
  }

  list(dat = dat, mr = mr_res, extras = extras, window_kb = window_kb, tag = tag)
}

primary <- run_window(500,  "0.5")
sens_1M <- run_window(1000, "1.0")

# ---- write ----
write_result <- function(r, suffix) {
  if (is.null(r)) return(invisible())
  write.csv(r$dat,
    file.path(OUT_DIR, sprintf("C_cis_TF_dat_%s.csv", suffix)), row.names=FALSE)
  write.csv(r$mr,
    file.path(OUT_DIR, sprintf("C_cis_TF_mr_%s.csv", suffix)), row.names=FALSE)
  if (!is.null(r$extras$egger_intercept))
    write.csv(r$extras$egger_intercept,
      file.path(OUT_DIR, sprintf("C_cis_TF_egger_%s.csv", suffix)), row.names=FALSE)
  if (!is.null(r$extras$Q))
    write.csv(r$extras$Q,
      file.path(OUT_DIR, sprintf("C_cis_TF_Q_%s.csv", suffix)), row.names=FALSE)
}
write_result(primary, "500kb")
write_result(sens_1M, "1Mb")

# ---- headline ----
summary_lines <- c("== Analysis C: Cis-only TF-region MR ==", "")
add_hd <- function(r, label) {
  if (is.null(r)) return(c(sprintf("%s: no result", label), ""))
  ivw <- r$mr[r$mr$method %in% c("Inverse variance weighted", "Wald ratio"), ]
  hdr <- c(sprintf("%s (n=%d SNPs):", label, nrow(r$dat)))
  for (i in seq_len(nrow(r$mr))) {
    hdr <- c(hdr, sprintf("  %-30s b=%.4f  se=%.4f  p=%.3e",
      r$mr$method[i], r$mr$b[i], r$mr$se[i], r$mr$pval[i]))
  }
  if (!is.null(r$extras$egger_intercept)) {
    hdr <- c(hdr, sprintf("  Egger intercept p=%.3f",
                          r$extras$egger_intercept$pval))
  }
  c(hdr, "")
}
summary_lines <- c(summary_lines,
                   add_hd(primary, "Primary: TF +/- 500 kb"),
                   add_hd(sens_1M, "Sensitivity: TF +/- 1 Mb"))
writeLines(summary_lines, file.path(OUT_DIR, "C_cis_TF_summary.txt"))
cat("\n", paste(summary_lines, collapse = "\n"), "\n", sep = "")
