##############################################################################
# R1_D_coloc_TF.R
#
# Reviewer response D: Bayesian colocalization at the TF locus between
# transferrin (ieu-a-1052) and circulating total-tau (ebi-a-GCST90095138).
#
# Pulls regional summary statistics for a chr3 window around the TF gene
# (default ±500 kb of TF gene body), then runs coloc::coloc.abf.
#
# H0: no signal in either
# H1: signal only in transferrin
# H2: signal only in tau
# H3: signals in both but distinct causal variants
# H4: signals in both driven by same causal variant  <-- what we want
##############################################################################

args <- commandArgs(trailingOnly = TRUE)
PROJECT_DIR <- if (length(args) >= 1) args[1] else normalizePath("..", mustWork = TRUE)
OUT_DIR     <- file.path(PROJECT_DIR, "results", "reviewer_response")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(ieugwasr); library(coloc); library(dplyr)
})
stopifnot(nchar(Sys.getenv("OPENGWAS_JWT")) > 100)

# TF gene (GRCh37): chr3:133,464,073 - 133,494,388  (Ensembl)
TF_CHR   <- 3L
TF_START <- 133464073L
TF_END   <- 133494388L
WINDOW_KB <- 500L
CHR_STR <- as.character(TF_CHR)
WIN_START <- TF_START - WINDOW_KB * 1000
WIN_END   <- TF_END   + WINDOW_KB * 1000

cat(sprintf("Coloc window: chr%s:%d-%d (%.2f Mb)\n",
            CHR_STR, WIN_START, WIN_END, (WIN_END - WIN_START) / 1e6))

# Sample sizes from source GWAS
N_TF  <- 23986    # Benyamin 2014 transferrin
N_TAU <- 14721    # Sarnowski 2022 total-tau

fetch_region <- function(id, chr, start, end, label) {
  cat(sprintf("  Fetching %s [%s] chr%s:%d-%d ...", label, id, chr, start, end))
  d <- tryCatch(
    associations(variants = sprintf("%s:%d-%d", chr, start, end),
                 id = id, proxies = 0),
    error = function(e) { cat(" ERR:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(d) || nrow(d) == 0) { cat(" no data\n"); return(NULL) }
  cat(sprintf(" %d SNPs\n", nrow(d)))
  # standardise
  d
}

tf  <- fetch_region("ieu-a-1052",         CHR_STR, WIN_START, WIN_END, "TF exposure")
tau <- fetch_region("ebi-a-GCST90095138", CHR_STR, WIN_START, WIN_END, "total-tau outcome")

if (is.null(tf) || is.null(tau)) {
  stop("Failed to fetch summary stats for coloc — cannot proceed.")
}

# Show schema in case coloc needs specific column names
cat("\ntransferrin columns:", paste(names(tf), collapse=", "), "\n")
cat("total-tau columns:  ", paste(names(tau), collapse=", "), "\n")

# The `associations` return usually has: rsid, chr, position, ea, nea, eaf, beta, se, p, n, ...
# Normalise
norm_df <- function(d) {
  need <- c("rsid","position","beta","se","p")
  ok <- need %in% names(d)
  if (!all(ok)) {
    stop("Missing expected columns: ",
         paste(need[!ok], collapse=", "),
         " — got: ", paste(names(d), collapse=", "))
  }
  eaf_num <- if ("eaf" %in% names(d)) as.numeric(d$eaf) else rep(NA_real_, nrow(d))
  data.frame(
    snp = d$rsid,
    position = as.integer(d$position),
    beta = as.numeric(d$beta),
    varbeta = as.numeric(d$se)^2,
    pvalue = as.numeric(d$p),
    MAF = pmin(eaf_num, 1 - eaf_num),
    stringsAsFactors = FALSE
  )
}
d_tf  <- norm_df(tf)
d_tau <- norm_df(tau)

# Benyamin 2014 transferrin (ieu-a-1052) does not deposit EAF. Borrow MAF from
# the outcome GWAS (both European, same SNPs).
if (all(is.na(d_tf$MAF)) && !all(is.na(d_tau$MAF))) {
  cat("Borrowing MAF from tau GWAS (transferrin GWAS has no EAF)\n")
  m <- match(d_tf$snp, d_tau$snp)
  d_tf$MAF <- d_tau$MAF[m]
}

# Drop invalid rows
d_tf  <- d_tf[is.finite(d_tf$beta) & is.finite(d_tf$varbeta) &
              d_tf$varbeta > 0 & is.finite(d_tf$MAF) & d_tf$MAF > 0.005, ]
d_tau <- d_tau[is.finite(d_tau$beta) & is.finite(d_tau$varbeta) &
               d_tau$varbeta > 0 & is.finite(d_tau$MAF) & d_tau$MAF > 0.005, ]

# Intersect on SNP id
common <- intersect(d_tf$snp, d_tau$snp)
cat(sprintf("\nSNPs after QC: TF %d, tau %d, common %d\n",
            nrow(d_tf), nrow(d_tau), length(common)))
if (length(common) < 20) {
  cat("WARNING: fewer than 20 common SNPs — coloc will be very noisy.\n")
}
d_tf  <- d_tf[match(common, d_tf$snp), ]
d_tau <- d_tau[match(common, d_tau$snp), ]

# ---- coloc.abf ----
D1 <- list(
  snp = d_tf$snp, position = d_tf$position, type = "quant",
  beta = d_tf$beta, varbeta = d_tf$varbeta,
  MAF = d_tf$MAF, N = N_TF
)
D2 <- list(
  snp = d_tau$snp, position = d_tau$position, type = "quant",
  beta = d_tau$beta, varbeta = d_tau$varbeta,
  MAF = d_tau$MAF, N = N_TAU
)

cat("\nRunning coloc.abf ...\n")
res <- coloc.abf(dataset1 = D1, dataset2 = D2,
                 p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
print(res$summary)

# Top posterior SNP
if (!is.null(res$results)) {
  top <- res$results[order(-res$results$SNP.PP.H4)[1], ]
  cat("\nTop SNP by PP.H4:\n"); print(top)
  write.csv(res$results,
            file.path(OUT_DIR, "D_coloc_TF_snp_posteriors.csv"),
            row.names = FALSE)
}

# Save summary
sum_df <- data.frame(
  window_kb = WINDOW_KB,
  n_common  = length(common),
  nsnps     = res$summary["nsnps"],
  PP.H0.abf = res$summary["PP.H0.abf"],
  PP.H1.abf = res$summary["PP.H1.abf"],
  PP.H2.abf = res$summary["PP.H2.abf"],
  PP.H3.abf = res$summary["PP.H3.abf"],
  PP.H4.abf = res$summary["PP.H4.abf"]
)
write.csv(sum_df, file.path(OUT_DIR, "D_coloc_TF_summary.csv"), row.names = FALSE)

lines <- c(
  "== Analysis D: coloc.abf at TF locus (transferrin ↔ total-tau) ==",
  sprintf("Window: chr%s:%d-%d (%.2f Mb)", CHR_STR, WIN_START, WIN_END,
          (WIN_END - WIN_START)/1e6),
  sprintf("Common SNPs after QC: %d", length(common)),
  sprintf("PP.H0 = %.3f  no signal", res$summary["PP.H0.abf"]),
  sprintf("PP.H1 = %.3f  transferrin only", res$summary["PP.H1.abf"]),
  sprintf("PP.H2 = %.3f  tau only",         res$summary["PP.H2.abf"]),
  sprintf("PP.H3 = %.3f  both, DIFFERENT causal variants", res$summary["PP.H3.abf"]),
  sprintf("PP.H4 = %.3f  both, SAME causal variant  <== colocalization", res$summary["PP.H4.abf"])
)
writeLines(lines, file.path(OUT_DIR, "D_coloc_TF_summary.txt"))
cat("\n", paste(lines, collapse="\n"), "\n", sep="")
