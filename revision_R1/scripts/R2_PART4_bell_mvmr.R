##############################################################################
# R2_PART4_bell_mvmr.R  (JAD revision round 2, Part 4 confirmatory)
#
# 2-way MVMR: transferrin (Benyamin ieu-a-1052) + serum iron (Bell 2021,
# GCST011367, n=163,511, LOCAL deCODE file) -> circulating total-tau (Sarnowski).
# Mirrors R1_E_MVMR_extended.R exactly; only the iron source changes from
# Benyamin (ieu-a-1049, overlapping) to Bell (independent samples).
# With independent exposure samples, gencov = 0 is EXACTLY justified.
#
# Prereq: extract ironserum.7z into data/bell2021_iron/ (e.g. via 7zr.exe).
# Auto-detects deCODE column names; prints the header first for verification.
#
# Writes: results/reviewer_response/TableS8_bell_mvmr.csv
#         results/reviewer_response/P4_bell_mvmr_summary.md
#         results/reviewer_response/P4_bell_mvmr_snps.csv
##############################################################################

readRenviron('C:/Users/ashwi/OneDrive/Documents/.Renviron')
.libPaths(c('C:/Users/ashwi/AppData/Local/R/win-library/4.6', .libPaths()))
suppressPackageStartupMessages({
  library(ieugwasr); library(TwoSampleMR); library(MVMR)
  have_dt <- requireNamespace("data.table", quietly = TRUE)
})
stopifnot(nchar(Sys.getenv("OPENGWAS_JWT")) > 100)

PROJECT_DIR <- "C:/Users/ashwi/OneDrive/Documents/EB1A Docs/MR Paper"
REV      <- file.path(PROJECT_DIR, "results", "reviewer_response")
BELL_DIR <- file.path(PROJECT_DIR, "data", "bell2021_iron")
TF_ID  <- "ieu-a-1052"; OUT_ID <- "ebi-a-GCST90095138"

# ---- locate + read the Bell iron sumstats ----
cand <- list.files(BELL_DIR, full.names = TRUE,
                   pattern = "\\.(txt|tsv|csv|gz|assoc|out)$", ignore.case = TRUE)
cand <- cand[order(-file.info(cand)$size)]
stopifnot("No Bell sumstats file found in data/bell2021_iron/" = length(cand) >= 1)
BELL_FILE <- cand[1]
cat("Bell file:", BELL_FILE, "\n")

bell <- if (have_dt) data.table::fread(BELL_FILE, showProgress = FALSE) else
        read.table(BELL_FILE, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
bell <- as.data.frame(bell)
cat("Bell rows:", nrow(bell), " cols:", paste(names(bell), collapse=", "), "\n")

# ---- auto-detect columns (deCODE / GWAS-Catalog variants) ----
pick <- function(cands) { m <- intersect(tolower(cands), tolower(names(bell)))
  if (!length(m)) return(NA_character_); names(bell)[match(m[1], tolower(names(bell)))] }
c_rsid <- pick(c("rsids","rsid","name","markername","snp","variant_id"))
c_ea   <- pick(c("effectallele","effect_allele","ea","alt","allele1","a1"))
c_oa   <- pick(c("otherallele","other_allele","oa","ref","allele0","allele2","a2"))
c_beta <- pick(c("beta","effect","effect_size","b"))
c_se   <- pick(c("se","standard_error","standarderror","sebeta"))
c_p    <- pick(c("pval","p","pvalue","p_value","p-value"))
c_eaf  <- pick(c("eaf","effect_allele_frequency","imputed_af","impmaf","maf","freq"))
c_n    <- pick(c("n","samplesize","sample_size"))
cat(sprintf("Detected: rsid=%s ea=%s oa=%s beta=%s se=%s p=%s eaf=%s n=%s\n",
            c_rsid,c_ea,c_oa,c_beta,c_se,c_p,c_eaf,c_n))
stopifnot("Could not detect required Bell columns" =
          !any(is.na(c(c_rsid,c_ea,c_oa,c_beta,c_se,c_p))))

bell_std <- data.frame(
  SNP  = bell[[c_rsid]], ea = toupper(bell[[c_ea]]), oa = toupper(bell[[c_oa]]),
  beta = as.numeric(bell[[c_beta]]), se = as.numeric(bell[[c_se]]),
  p    = as.numeric(bell[[c_p]]),
  eaf  = if (!is.na(c_eaf)) as.numeric(bell[[c_eaf]]) else NA,
  n    = if (!is.na(c_n))   as.numeric(bell[[c_n]])   else 163511,
  stringsAsFactors = FALSE)
# some deCODE files list multiple rsids per row ("rs1,rs2"); take first
bell_std$SNP <- sub(",.*$", "", bell_std$SNP)

# ---- Bell iron instruments: p<5e-8, clump r2<0.001 (EUR) ----
sig <- bell_std[!is.na(bell_std$p) & bell_std$p < 5e-8 & grepl("^rs", bell_std$SNP), ]
cat("Bell iron p<5e-8 SNPs:", nrow(sig), "\n")
fe_clump <- ld_clump(dplyr::tibble(rsid = sig$SNP, pval = sig$p),
                     clump_r2 = 0.001, clump_kb = 10000, pop = "EUR",
                     opengwas_jwt = Sys.getenv("OPENGWAS_JWT"))
fe_inst <- fe_clump$rsid
cat("Bell iron independent instruments:", length(fe_inst), "\n")

# ---- transferrin instruments (unchanged from primary) ----
tf_inst <- extract_instruments(TF_ID, p1 = 5e-8, r2 = 0.001, kb = 10000)$SNP
union_snps <- unique(c(tf_inst, fe_inst))
cat("Union:", length(union_snps), "SNPs\n")

# ---- SNP effects on each node ----
fetch <- function(snps, id) {
  a <- associations(variants = snps, id = id, proxies = 0)
  data.frame(SNP=a$rsid, beta=a$beta, se=a$se, ea=a$ea, oa=a$nea,
             p=a$p, stringsAsFactors=FALSE)
}
tf_e  <- fetch(union_snps, TF_ID)
tau_e <- fetch(union_snps, OUT_ID)
fe_e  <- bell_std[match(union_snps, bell_std$SNP), c("SNP","beta","se","ea","oa","p")]
fe_e  <- fe_e[!is.na(fe_e$SNP), ]

common <- Reduce(intersect, list(tf_e$SNP, fe_e$SNP, tau_e$SNP))
cat("Common across transferrin + Bell iron + tau:", length(common), "SNPs\n")
g <- function(df,s) df[match(s,df$SNP),]
tf_r <- g(tf_e,common); fe_r <- g(fe_e,common); tau_r <- g(tau_e,common)

# ---- align iron + tau to transferrin effect allele ----
al <- function(ea,oa,b,rea,roa){ flip<-(ea==roa)&(oa==rea); same<-(ea==rea)&(oa==roa)
  b[flip]<--b[flip]; list(b=b, keep=flip|same) }
afe  <- al(fe_r$ea, fe_r$oa, fe_r$beta, tf_r$ea, tf_r$oa)
atau <- al(tau_r$ea,tau_r$oa,tau_r$beta,tf_r$ea, tf_r$oa)
keep <- afe$keep & atau$keep
cat("Allele-alignable:", sum(keep), "/", length(common), "\n")
tf_r<-tf_r[keep,]; fe_r<-fe_r[keep,]; tau_r<-tau_r[keep,]
common<-common[keep]; afe$b<-afe$b[keep]; atau$b<-atau$b[keep]

# ---- MVMR (gencov=0 EXACT here: independent exposure samples) ----
bx <- cbind(tf_r$beta, afe$b); bxse <- cbind(tf_r$se, fe_r$se)
colnames(bx)<-colnames(bxse)<-c("transferrin","serum_iron_bell")
mvin  <- format_mvmr(BXGs=bx, seBXGs=bxse, BYG=atau$b, seBYG=tau_r$se, RSID=common)
Fc    <- strength_mvmr(r_input=mvin, gencov=0)
res   <- ivw_mvmr(r_input=mvin)
QA    <- pleiotropy_mvmr(r_input=mvin, gencov=0)
r_exp <- cor(bx[,1], bx[,2])

out <- data.frame(exposure=c("transferrin","serum_iron_bell"),
  beta=res[,"Estimate"], se=res[,"Std. Error"], t=res[,"t value"],
  pval=res[,"Pr(>|t|)"], cond_F=as.numeric(Fc[1,]))
write.csv(out, file.path(REV,"TableS8_bell_mvmr.csv"), row.names=FALSE)
write.csv(data.frame(SNP=common, beta_tf=bx[,1], se_tf=bxse[,1],
  beta_fe=bx[,2], se_fe=bxse[,2], beta_tau=atau$b, se_tau=tau_r$se),
  file.path(REV,"P4_bell_mvmr_snps.csv"), row.names=FALSE)

md <- c("# Part 4 confirmatory — Bell 2021 iron MVMR (independent exposure samples)","",
  sprintf("Transferrin (Benyamin ieu-a-1052) + serum iron (Bell 2021 GCST011367, n=163,511, deCODE local) -> total-tau."),
  sprintf("Instruments: %d transferrin + %d Bell-iron -> union, %d after alignment+intersection. gencov=0 EXACT (independent samples; exposure-exposure overlap = 0).", length(tf_inst), length(fe_inst), nrow(bx)),"",
  sprintf("- Transferrin: b=%.4f, SE=%.4f, p=%.3g (t-based, df=%d), conditional F=%.1f", out$beta[1],out$se[1],out$pval[1], nrow(bx)-2, out$cond_F[1]),
  sprintf("- Serum iron (Bell): b=%.4f, SE=%.4f, p=%.3g, conditional F=%.1f", out$beta[2],out$se[2],out$pval[2], out$cond_F[2]),
  sprintf("- MVMR Q_A = %.2f (p=%.3g); SNP-effect r(transferrin,iron) = %.3f", as.numeric(QA$Qstat), as.numeric(QA$Qpval), r_exp),"",
  "Compare to the Benyamin-iron MVMR (Analysis E, overlapping samples): transferrin b=-0.0437, p=1.4e-3, cond-F=278.",
  "If transferrin holds here with an INDEPENDENT iron GWAS, the finding is robust to the exposure-overlap concern.")
writeLines(md, file.path(REV,"P4_bell_mvmr_summary.md"))
cat("\n", paste(md, collapse="\n"), "\n", sep="")
