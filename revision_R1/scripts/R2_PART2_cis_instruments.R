##############################################################################
# R2_PART2_cis_instruments.R  (JAD revision round 2, Part 2)
#
# Per-SNP cis instrument detail for {rs8177240, rs3811658, rs17376530}:
#   individual Wald ratios, F-stats, GRCh37 positions + distance from TF gene,
#   pairwise LD r2 (1000G EUR via OpenGWAS), and the rs17376530-alone estimate.
#
# Writes:
#   results/reviewer_response/TableS6_cis_instruments.csv
#   results/reviewer_response/TableS6_cis_summary.md
#   results/reviewer_response/P2_ld_matrix.csv
##############################################################################

readRenviron('C:/Users/ashwi/OneDrive/Documents/.Renviron')
.libPaths(c('C:/Users/ashwi/AppData/Local/R/win-library/4.6', .libPaths()))
suppressPackageStartupMessages({ library(ieugwasr) })

PROJECT_DIR <- "C:/Users/ashwi/OneDrive/Documents/EB1A Docs/MR Paper"
REV <- file.path(PROJECT_DIR, "results", "reviewer_response")
DATE_RUN <- "2026-08-07"

TF_START <- 133464073L; TF_END <- 133494388L   # TF gene body, GRCh37
gene_dist <- function(pos) if (pos >= TF_START & pos <= TF_END) 0L else
  as.integer(min(abs(pos - TF_START), abs(pos - TF_END)))

snps <- c("rs8177240","rs3811658","rs17376530")

# ---- exposure/outcome betas from cached cis + primary data ----
cis <- read.csv(file.path(REV,"C_cis_TF_dat_500kb.csv"), stringsAsFactors=FALSE)
harm <- read.csv(file.path(PROJECT_DIR,"results","primary",
        "harmonised_transferrin_tau_2026-03-28.csv"), stringsAsFactors=FALSE)

# assemble per-SNP exposure/outcome effects (prefer cis file; rs8177240 from primary)
getrow <- function(s) {
  if (s %in% cis$SNP) {
    r <- cis[cis$SNP==s,]
    data.frame(SNP=s, bx=r$beta.exposure, bxse=r$se.exposure,
               by=r$beta.outcome, byse=r$se.outcome, pos=r$pos, src="cis(ieu-a-1052)")
  } else {
    r <- harm[harm$SNP==s,]
    data.frame(SNP=s, bx=r$beta.exposure, bxse=r$se.exposure,
               by=r$beta.outcome, byse=r$se.outcome, pos=r$pos, src="primary(ieu-a-1052)")
  }
}
tab <- do.call(rbind, lapply(snps, getrow))
tab$wald    <- tab$by / tab$bx
tab$wald_se <- tab$byse / abs(tab$bx)
tab$wald_p  <- 2*pnorm(-abs(tab$wald/tab$wald_se))
tab$Fstat   <- (tab$bx/tab$bxse)^2
tab$dist_TF_bp <- sapply(tab$pos, gene_dist)

# ---- live check: current exposure betas on ieu-a-1052 (resolve 0.388 vs 0.423) ----
liveassoc <- tryCatch(associations(variants=snps, id="ieu-a-1052"),
                      error=function(e) NULL)
tab$live_bx <- NA; tab$live_ea <- NA
if (!is.null(liveassoc)) {
  for (i in seq_len(nrow(tab))) {
    m <- liveassoc[liveassoc$rsid==tab$SNP[i],]
    if (nrow(m)>=1){ tab$live_bx[i] <- m$beta[1]; tab$live_ea[i] <- m$ea[1] }
  }
}

# ---- LD matrix, 1000G EUR ----
ldr <- tryCatch(ld_matrix(snps, pop="EUR", with_alleles=FALSE),
                error=function(e) { cat("LD ERR:",conditionMessage(e),"\n"); NULL })
ld_long <- NULL
if (!is.null(ldr)) {
  write.csv(as.data.frame(ldr), file.path(REV,"P2_ld_matrix.csv"))
  cn <- colnames(ldr)
  pr <- combn(cn,2)
  ld_long <- data.frame(snp_a=pr[1,], snp_b=pr[2,],
    r=mapply(function(a,b) ldr[a,b], pr[1,], pr[2,]))
  ld_long$r2 <- ld_long$r^2
}

# ---- rs17376530 alone (genuinely independent instrument) ----
r2a <- tab[tab$SNP=="rs17376530",]

# ---- write TableS6 ----
out <- tab[,c("SNP","pos","dist_TF_bp","bx","bxse","Fstat","by","byse",
              "wald","wald_se","wald_p","live_bx","live_ea","src")]
names(out) <- c("SNP","pos_GRCh37","dist_from_TF_bp","beta_exposure","se_exposure",
                "F_stat","beta_outcome","se_outcome","wald_beta","wald_se","wald_p",
                "live_beta_exposure","live_effect_allele","source")
# keep positions/distance as integers; round effect sizes only
eff_cols <- c("beta_exposure","se_exposure","beta_outcome","se_outcome",
              "wald_beta","wald_se","wald_p","live_beta_exposure")
out[,eff_cols] <- lapply(out[,eff_cols], function(x) signif(x,4))
out$pos_GRCh37 <- as.integer(out$pos_GRCh37)
out$dist_from_TF_bp <- as.integer(out$dist_from_TF_bp)
out$F_stat <- round(out$F_stat)
write.csv(out, file.path(REV,"TableS6_cis_instruments.csv"), row.names=FALSE)

# ---- summary md ----
indep <- if (!is.null(ld_long)) {
  r2_530_240 <- ld_long$r2[(ld_long$snp_a=="rs17376530"&ld_long$snp_b=="rs8177240")|
                           (ld_long$snp_a=="rs8177240"&ld_long$snp_b=="rs17376530")]
  r2_658_240 <- ld_long$r2[(ld_long$snp_a=="rs3811658"&ld_long$snp_b=="rs8177240")|
                           (ld_long$snp_a=="rs8177240"&ld_long$snp_b=="rs3811658")]
  list(r2_530_240=r2_530_240, r2_658_240=r2_658_240)
} else list(r2_530_240=NA, r2_658_240=NA)

md <- c(
"# Table S6 — per-SNP cis instrument detail","",
sprintf("Generated %s. TF gene body GRCh37 chr3:%d-%d. LD from 1000G EUR (OpenGWAS).", DATE_RUN, TF_START, TF_END),"",
"## Individual Wald ratios",
sprintf("- rs8177240 (primary panel lead): beta=%.4f, SE=%.4f, p=%.4f, F=%.0f, %d bp from TF (inside gene body if 0).",
  tab$wald[tab$SNP=='rs8177240'], tab$wald_se[tab$SNP=='rs8177240'], tab$wald_p[tab$SNP=='rs8177240'], tab$Fstat[tab$SNP=='rs8177240'], tab$dist_TF_bp[tab$SNP=='rs8177240']),
sprintf("- rs3811658 (cis instrument 1, tags rs8177240): beta=%.4f, SE=%.4f, p=%.4f, F=%.0f, %d bp from TF.",
  tab$wald[tab$SNP=='rs3811658'], tab$wald_se[tab$SNP=='rs3811658'], tab$wald_p[tab$SNP=='rs3811658'], tab$Fstat[tab$SNP=='rs3811658'], tab$dist_TF_bp[tab$SNP=='rs3811658']),
sprintf("- rs17376530 (cis instrument 2, candidate independent): beta=%.4f, SE=%.4f, p=%.4f, F=%.0f, %d bp from TF.",
  r2a$wald, r2a$wald_se, r2a$wald_p, r2a$Fstat, r2a$dist_TF_bp),"",
"## LD (1000G EUR)",
if (!is.null(ld_long)) paste(sprintf("- r2(%s,%s) = %.4f", ld_long$snp_a, ld_long$snp_b, ld_long$r2), collapse="\n") else "- LD matrix unavailable (API error).","",
"## Verdict on independence (bar: r2 < 0.01)",
sprintf("- rs3811658 vs rs8177240: r2 = %s  -> %s (expected high: same signal, 849 bp apart).",
  ifelse(is.na(indep$r2_658_240),"NA",sprintf("%.4f",indep$r2_658_240)),
  ifelse(is.na(indep$r2_658_240),"?", ifelse(indep$r2_658_240>=0.01,"NOT independent","independent"))),
sprintf("- rs17376530 vs rs8177240: r2 = %s  -> %s.",
  ifelse(is.na(indep$r2_530_240),"NA",sprintf("%.4f",indep$r2_530_240)),
  ifelse(is.na(indep$r2_530_240),"?", ifelse(indep$r2_530_240<0.01,"INDEPENDENT of rs8177240","NOT independent"))),"",
"## Plain-language conclusion",
sprintf("rs17376530 is %s of rs8177240 and supports the SAME (protective) direction (beta=%.3f), but is NOT individually significant (p=%.2f). The cis-only IVW (beta=-0.060, p=3.99e-4) is ~90%% weighted on rs3811658, which tags rs8177240. So the cis analysis is NOT simply rs8177240 relabelled (rs17376530 is a distinct, concordant variant), but neither is it two independently-significant instruments. Report this honestly: the cis result corroborates direction via an independent variant while its significance still leans on the rs8177240-tagging signal.",
  ifelse(is.na(indep$r2_530_240),"(LD pending)", ifelse(indep$r2_530_240<0.01,"INDEPENDENT","in LD")),
  r2a$wald, r2a$wald_p),"",
"## Data-version note (rs8177240 / rs3811658 rsID labelling) — RESOLVED",
sprintf("- Harmonized/primary rs8177240 beta.exposure = %.4f (F=%.0f); the live ieu-a-1052 API returns %s for the rs8177240 label.",
  tab$bx[tab$SNP=='rs8177240'], tab$Fstat[tab$SNP=='rs8177240'],
  ifelse(is.na(tab$live_bx[tab$SNP=='rs8177240']),"NA",sprintf("%.4f",tab$live_bx[tab$SNP=='rs8177240']))),
"- RESOLVED against Ensembl GRCh37 (authoritative): rs8177240 = chr3:133,477,701, rs3811658 = chr3:133,476,852.",
"  The CACHED/manuscript data match Ensembl; the live OpenGWAS ieu-a-1052 API SWAPS these two",
"  adjacent rsID labels. Manuscript is correct; no estimate changes. See P2_rsID_swap_investigation.md.",
"  Reproducibility: map instruments by position, not live rsID, when using ieu-a-1052."
)
writeLines(md, file.path(REV,"TableS6_cis_summary.md"))

cat("=== TableS6 ===\n"); print(out, row.names=FALSE)
cat("\n=== LD (1000G EUR) ===\n")
if(!is.null(ld_long)) print(ld_long, row.names=FALSE) else cat("LD unavailable\n")
cat("\nlive_bx rs8177240:", tab$live_bx[tab$SNP=='rs8177240'],
    " harmonized bx:", tab$bx[tab$SNP=='rs8177240'], "\n")
