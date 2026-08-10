##############################################################################
# R1_B_leave_locus_out.R
#
# Reviewer response B: Leave-locus-out sensitivity for the primary transferrin
# -> circulating total-tau IVW estimate.
#
# Loci to drop: HLA (rs9268633), FADS (rs174577), HFE (rs1800562),
#               NAT2 (rs1495741), plus all pairwise + joint combinations.
#
# Also reports base-R IVW, weighted median, and MR-Egger for each subset,
# and Cochran's Q heterogeneity.
#
# Reads:  results/primary/harmonised_transferrin_tau_2026-03-28.csv
# Writes: results/reviewer_response/B_leave_locus_out.csv
#         results/reviewer_response/B_leave_locus_out_summary.txt
#
# No API, no external packages beyond base R.
##############################################################################

args <- commandArgs(trailingOnly = TRUE)
PROJECT_DIR <- if (length(args) >= 1) args[1] else normalizePath("..", mustWork = TRUE)
OUT_DIR     <- file.path(PROJECT_DIR, "results", "reviewer_response")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

dat <- read.csv(file.path(PROJECT_DIR, "results", "primary",
                          "harmonised_transferrin_tau_2026-03-28.csv"),
                stringsAsFactors = FALSE)
dat <- dat[dat$mr_keep == TRUE, ]
dat$wald    <- dat$beta.outcome / dat$beta.exposure
dat$wald_se <- dat$se.outcome / abs(dat$beta.exposure)

# Locus map from manuscript / SNP annotations
locus_map <- list(
  TF   = c("rs8177240", "rs7646473", "rs9990333"),      # chr3 transferrin region + downstream
  HFE  = c("rs1800562"),                                 # chr6 HFE
  HLA  = c("rs9268633"),                                 # chr6 HLA
  FADS = c("rs174577"),                                  # chr11 FADS1/2
  NAT2 = c("rs1495741"),                                 # chr8 NAT2
  TFRC = c("rs744653")                                   # chr3 TFRC region (per annotation)
)

# --- helpers ---
ivw_fixed <- function(b, s) {
  w  <- 1 / s^2
  bh <- sum(w * b) / sum(w)
  sh <- sqrt(1 / sum(w))
  z  <- bh / sh
  p  <- 2 * pnorm(-abs(z))
  Q  <- sum(w * (b - bh)^2)
  Qp <- pchisq(Q, df = length(b) - 1, lower.tail = FALSE)
  I2 <- max(0, (Q - (length(b) - 1)) / Q * 100)
  list(beta = bh, se = sh, pval = p, Q = Q, Q_p = Qp, I2 = I2, k = length(b))
}

weighted_median <- function(b, s) {
  # Bowden 2016 weighted median with inverse-variance weights.
  w  <- 1 / s^2
  ord <- order(b)
  b_s <- b[ord]; w_s <- w[ord]
  cw <- cumsum(w_s) - 0.5 * w_s
  cw <- cw / sum(w_s)
  # locate median
  hi <- which(cw >= 0.5)[1]
  lo <- hi - 1
  if (lo < 1) return(list(beta = b_s[hi], se = NA_real_))
  bh <- b_s[lo] + (b_s[hi] - b_s[lo]) * (0.5 - cw[lo]) / (cw[hi] - cw[lo])
  # bootstrap SE
  set.seed(1); B <- 1000
  boots <- replicate(B, {
    bb <- rnorm(length(b), b, s)
    ord <- order(bb)
    bb <- bb[ord]; ww <- w[ord]
    cw <- cumsum(ww) - 0.5 * ww
    cw <- cw / sum(ww)
    hi <- which(cw >= 0.5)[1]; lo <- hi - 1
    if (lo < 1) bb[hi] else
      bb[lo] + (bb[hi] - bb[lo]) * (0.5 - cw[lo]) / (cw[hi] - cw[lo])
  })
  sh <- sd(boots)
  z <- bh / sh; p <- 2 * pnorm(-abs(z))
  list(beta = bh, se = sh, pval = p)
}

mr_egger <- function(b, s, bx) {
  # Standard MR-Egger: regress Wald ratios on 1/|bx| weighted by 1/s^2.
  # Equivalent to weighted regression of beta.outcome on beta.exposure with intercept.
  # Uses the harmonized inputs (aligns effect allele so bx is positive per TwoSampleMR harmonise action=2).
  # We regress beta.outcome ~ beta.exposure weighted by 1/se.outcome^2.
  if (length(b) < 3) return(list(beta = NA, se = NA, pval = NA,
                                  intercept = NA, int_se = NA, int_p = NA))
  # b here is Wald ratio; convert back with bx
  by <- b * bx
  bys <- s * abs(bx)
  w <- 1 / bys^2
  fit <- lm(by ~ bx, weights = w)
  co  <- summary(fit)$coefficients
  list(
    beta = co["bx", "Estimate"], se = co["bx", "Std. Error"],
    pval = co["bx", "Pr(>|t|)"],
    intercept = co["(Intercept)", "Estimate"],
    int_se = co["(Intercept)", "Std. Error"],
    int_p  = co["(Intercept)", "Pr(>|t|)"]
  )
}

run_subset <- function(df, label) {
  bx <- df$beta.exposure
  b  <- df$wald
  s  <- df$wald_se
  iv <- ivw_fixed(b, s)
  wm <- if (length(b) >= 3) weighted_median(b, s) else list(beta=NA, se=NA, pval=NA)
  eg <- if (length(b) >= 3) mr_egger(b, s, bx) else
    list(beta=NA, se=NA, pval=NA, intercept=NA, int_se=NA, int_p=NA)
  data.frame(
    subset = label, n_snps = iv$k,
    ivw_beta = iv$beta, ivw_se = iv$se, ivw_pval = iv$pval,
    ivw_Q = iv$Q, ivw_Q_p = iv$Q_p, ivw_I2 = iv$I2,
    wm_beta  = wm$beta,  wm_se  = wm$se,  wm_pval  = wm$pval,
    egg_beta = eg$beta, egg_se = eg$se, egg_pval = eg$pval,
    egg_intercept = eg$intercept, egg_int_se = eg$int_se, egg_int_p = eg$int_p
  )
}

# --- run all subsets ---
subsets <- list(
  "all_8"                               = dat$SNP,
  "excl_rs8177240"                      = setdiff(dat$SNP, "rs8177240"),
  "excl_HLA"                            = setdiff(dat$SNP, locus_map$HLA),
  "excl_FADS"                           = setdiff(dat$SNP, locus_map$FADS),
  "excl_HFE"                            = setdiff(dat$SNP, locus_map$HFE),
  "excl_NAT2"                           = setdiff(dat$SNP, locus_map$NAT2),
  "excl_HLA_FADS"                       = setdiff(dat$SNP, c(locus_map$HLA, locus_map$FADS)),
  "excl_HLA_FADS_HFE_NAT2 (joint)"      = setdiff(dat$SNP, c(locus_map$HLA, locus_map$FADS,
                                                              locus_map$HFE, locus_map$NAT2)),
  "TF_only (rs8177240)"                 = "rs8177240",
  "cis_TF_3snps (rs8177240,rs9990333,rs7646473)" =
    c("rs8177240", "rs9990333", "rs7646473")
)

res <- do.call(rbind, lapply(names(subsets), function(nm) {
  keep <- dat[dat$SNP %in% subsets[[nm]], ]
  if (nrow(keep) == 0) return(NULL)
  run_subset(keep, nm)
}))

# ---- format ----
num_cols <- sapply(res, is.numeric)
res_out <- res
res_out[, num_cols] <- lapply(res_out[, num_cols], function(x) signif(x, 4))

write.csv(res_out, file.path(OUT_DIR, "B_leave_locus_out.csv"), row.names = FALSE)

# ---- short summary ----
summary_lines <- c(
  "== Analysis B: Leave-locus-out sensitivity ==",
  "",
  "Locus map:",
  "  TF   locus  -> rs8177240, rs7646473, rs9990333 (chr3)",
  "  HFE         -> rs1800562  (chr6)",
  "  HLA         -> rs9268633  (chr6)",
  "  FADS        -> rs174577   (chr11)",
  "  NAT2        -> rs1495741  (chr8)",
  "  TFRC region -> rs744653   (chr2, per SNP annotation table)",
  ""
)
tab <- res[, c("subset", "n_snps", "ivw_beta", "ivw_se", "ivw_pval",
               "wm_beta", "wm_pval", "egg_beta", "egg_pval", "egg_int_p",
               "ivw_Q_p", "ivw_I2")]
tab[, sapply(tab, is.numeric)] <- lapply(
  tab[, sapply(tab, is.numeric)], function(x) signif(x, 3))
summary_lines <- c(summary_lines, capture.output(print(tab, row.names = FALSE)))

writeLines(summary_lines, file.path(OUT_DIR, "B_leave_locus_out_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")
