##############################################################################
# R2_A4_window_counts.R
#
# Completes the one cell of Analysis A4 that could not be recovered offline:
# the number of variants in the colocalization window for each trait BEFORE
# harmonisation. Replicates the fetch in R1_D_coloc_TF.R exactly.
#
# Success criterion: the harmonised intersection must reproduce 696, the count
# used in the published colocalization. If it does not, the fetch does not
# match the R1 analysis and the counts must not be reported.
##############################################################################

.libPaths(c("C:/Users/ashwi/AppData/Local/R/win-library/4.6", .libPaths()))
suppressPackageStartupMessages({ library(ieugwasr) })

PROJECT_DIR <- "C:/Users/ashwi/OneDrive/Documents/EB1A Docs/MR Paper"
OUT_DIR     <- file.path(PROJECT_DIR, "results", "R2")

stopifnot(nchar(Sys.getenv("OPENGWAS_JWT")) > 100)

# identical to R1_D_coloc_TF.R
TF_START <- 133464073L; TF_END <- 133494388L; WINDOW_KB <- 500L
WIN_START <- TF_START - WINDOW_KB * 1000
WIN_END   <- TF_END   + WINDOW_KB * 1000
REGION <- sprintf("3:%d-%d", WIN_START, WIN_END)
cat(sprintf("Window: chr%s  (%.2f Mb)\n\n", REGION, (WIN_END - WIN_START) / 1e6))

fetch <- function(id, label) {
  cat(sprintf("Fetching %-12s [%s] ... ", label, id))
  d <- tryCatch(associations(variants = REGION, id = id, proxies = 0),
                error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL })
  if (is.null(d)) return(NULL)
  cat(sprintf("%d variants\n", nrow(d)))
  d
}

tf  <- fetch("ieu-a-1052",         "transferrin")
tau <- fetch("ebi-a-GCST90095138", "total-tau")
stopifnot(!is.null(tf), !is.null(tau))

# same QC as R1
norm <- function(d) data.frame(
  snp = d$rsid, beta = as.numeric(d$beta), varbeta = as.numeric(d$se)^2,
  MAF = { e <- if ("eaf" %in% names(d)) as.numeric(d$eaf) else NA_real_
          pmin(e, 1 - e) }, stringsAsFactors = FALSE)
d_tf <- norm(tf); d_tau <- norm(tau)
if (all(is.na(d_tf$MAF))) d_tf$MAF <- d_tau$MAF[match(d_tf$snp, d_tau$snp)]
ok <- function(d) d[is.finite(d$beta) & is.finite(d$varbeta) & d$varbeta > 0 &
                    is.finite(d$MAF) & d$MAF > 0.005, ]
q_tf <- ok(d_tf); q_tau <- ok(d_tau)
common <- intersect(q_tf$snp, q_tau$snp)

cat("\n--- counts ---\n")
cat(sprintf("  transferrin, in window, before harmonisation : %d\n", nrow(tf)))
cat(sprintf("  total-tau,   in window, before harmonisation : %d\n", nrow(tau)))
cat(sprintf("  transferrin, after QC (MAF>0.005, finite SE) : %d\n", nrow(q_tf)))
cat(sprintf("  total-tau,   after QC                        : %d\n", nrow(q_tau)))
cat(sprintf("  shared and analysed                          : %d\n", length(common)))
cat(sprintf("\n  published value was 696 -> %s\n",
            ifelse(length(common) == 696, "REPRODUCED",
                   "*** MISMATCH - DO NOT REPORT THESE COUNTS ***")))

res <- data.frame(
  item = c("Variants in window, transferrin GWAS, before harmonisation",
           "Variants in window, total-tau GWAS, before harmonisation",
           "Variants in window, transferrin GWAS, after QC",
           "Variants in window, total-tau GWAS, after QC",
           "Variants shared and analysed after harmonisation and QC"),
  value = c(nrow(tf), nrow(tau), nrow(q_tf), nrow(q_tau), length(common)),
  stringsAsFactors = FALSE)
p <- file.path(OUT_DIR, "coloc_window_counts.csv")
write.csv(res, p, row.names = FALSE)
cat(sprintf("\n  WROTE md5 %s  %s\n", tools::md5sum(p), basename(p)))
print(read.csv(p), row.names = FALSE)
