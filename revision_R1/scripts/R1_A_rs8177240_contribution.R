##############################################################################
# R1_A_rs8177240_contribution.R
#
# Reviewer response A: Explicit Wald ratio and IVW weight contribution of
# rs8177240 for the transferrin -> circulating total-tau primary analysis.
#
# Reads:  results/primary/harmonised_transferrin_tau_2026-03-28.csv
# Writes: results/reviewer_response/A_rs8177240_contribution.csv
#         results/reviewer_response/A_rs8177240_summary.txt
#
# No API, no external packages beyond base R.
##############################################################################

args <- commandArgs(trailingOnly = TRUE)
PROJECT_DIR <- if (length(args) >= 1) args[1] else
  normalizePath("..", mustWork = TRUE)
OUT_DIR     <- file.path(PROJECT_DIR, "results", "reviewer_response")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

dat <- read.csv(file.path(PROJECT_DIR, "results", "primary",
                          "harmonised_transferrin_tau_2026-03-28.csv"),
                stringsAsFactors = FALSE)
dat <- dat[dat$mr_keep == TRUE, ]
stopifnot(nrow(dat) == 8)

# ---- Per-SNP Wald ratio, SE, IVW weight ----
dat$wald   <- dat$beta.outcome / dat$beta.exposure
# Standard MR summary-stats Wald SE (delta-method, ignoring exposure SE per Burgess)
dat$wald_se <- dat$se.outcome / abs(dat$beta.exposure)
dat$wald_z  <- dat$wald / dat$wald_se
dat$wald_p  <- 2 * pnorm(-abs(dat$wald_z))
# IVW weight (inverse variance of Wald ratio)
dat$weight  <- 1 / dat$wald_se^2
dat$weight_pct <- 100 * dat$weight / sum(dat$weight)

per_snp <- dat[order(-dat$weight_pct),
               c("SNP", "chr", "pos", "beta.exposure", "se.exposure",
                 "beta.outcome", "se.outcome", "F_stat",
                 "wald", "wald_se", "wald_p", "weight", "weight_pct")]
per_snp[, sapply(per_snp, is.numeric)] <-
  lapply(per_snp[, sapply(per_snp, is.numeric)],
         function(x) signif(x, 4))

# ---- IVW fixed effects, full 8 SNPs ----
ivw_fixed <- function(beta, se) {
  w <- 1 / se^2
  b <- sum(w * beta) / sum(w)
  s <- sqrt(1 / sum(w))
  z <- b / s
  p <- 2 * pnorm(-abs(z))
  list(beta = b, se = s, z = z, pval = p, k = length(beta))
}

full <- ivw_fixed(dat$wald, dat$wald_se)
sub  <- ivw_fixed(dat$wald[dat$SNP != "rs8177240"],
                  dat$wald_se[dat$SNP != "rs8177240"])
wald_alone <- dat[dat$SNP == "rs8177240",
                  c("wald", "wald_se", "wald_p", "weight_pct")]

# ---- Write outputs ----
write.csv(per_snp, file.path(OUT_DIR, "A_rs8177240_contribution.csv"),
          row.names = FALSE)

summary_lines <- c(
  "== Analysis A: rs8177240 contribution to primary IVW estimate ==",
  "",
  sprintf("Primary IVW (all 8 SNPs):       beta = %.4f, SE = %.4f, p = %.3e",
          full$beta, full$se, full$pval),
  sprintf("IVW excluding rs8177240 (n=7):  beta = %.4f, SE = %.4f, p = %.3e",
          sub$beta, sub$se, sub$pval),
  "",
  "-- rs8177240 (TF gene) single-instrument statistics --",
  sprintf("Wald ratio            : %.4f (SE %.4f, p = %.3e)",
          wald_alone$wald, wald_alone$wald_se, wald_alone$wald_p),
  sprintf("F-statistic           : %.1f",
          dat$F_stat[dat$SNP == "rs8177240"]),
  sprintf("Share of IVW weight   : %.1f%%", wald_alone$weight_pct),
  "",
  "Per-SNP breakdown (ranked by IVW weight%):",
  ""
)
tbl <- capture.output(print(per_snp[, c("SNP", "chr", "F_stat",
                                        "wald", "wald_se", "wald_p",
                                        "weight_pct")],
                             row.names = FALSE))
summary_lines <- c(summary_lines, tbl)

writeLines(summary_lines,
           file.path(OUT_DIR, "A_rs8177240_summary.txt"))

cat(paste(summary_lines, collapse = "\n"), "\n")
