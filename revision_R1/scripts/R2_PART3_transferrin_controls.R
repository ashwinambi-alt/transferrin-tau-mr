##############################################################################
# R2_PART3_transferrin_controls.R  (JAD revision round 2, Part 3)
#
# Transferrin-specific positive/negative controls + DAG (mediator) check,
# using the SAME 8-SNP transferrin instrument set and the SAME TwoSampleMR
# pipeline as the primary analysis. Reviewer 2 noted the existing serum-iron
# controls validate the pipeline but NOT the transferrin instruments.
#
# Positive control : transferrin -> transferrin saturation (expect NEGATIVE)
# Negative controls: transferrin -> height ; transferrin -> educational attainment (expect null)
# DAG / mediator   : transferrin -> serum iron  (+ Steiger direction)
#
# NB: TIBC positive control is NOT in OpenGWAS (Bell TIBC = deCODE only); noted, skipped.
#
# Writes: results/reviewer_response/TableS7_transferrin_controls.csv
#         results/reviewer_response/TableS7_transferrin_controls_summary.md
##############################################################################

readRenviron('C:/Users/ashwi/OneDrive/Documents/.Renviron')
.libPaths(c('C:/Users/ashwi/AppData/Local/R/win-library/4.6', .libPaths()))
suppressPackageStartupMessages({ library(TwoSampleMR) })

PROJECT_DIR <- "C:/Users/ashwi/OneDrive/Documents/EB1A Docs/MR Paper"
REV <- file.path(PROJECT_DIR, "results", "reviewer_response")
DATE_RUN <- "2026-08-08"

# ---- Exposure: 8-SNP transferrin instrument set from the CACHED (Ensembl-consistent)
#      harmonized primary. Using cached effects avoids the live-API rs8177240/rs3811658
#      label swap (see P2_rsID_swap_investigation.md). ----
harm <- read.csv(file.path(PROJECT_DIR,"results","primary",
        "harmonised_transferrin_tau_2026-03-28.csv"), stringsAsFactors=FALSE)
harm <- harm[harm$mr_keep==TRUE,]
exp_raw <- data.frame(
  SNP           = harm$SNP,
  beta          = harm$beta.exposure,
  se            = harm$se.exposure,
  effect_allele = harm$effect_allele.exposure,
  other_allele  = harm$other_allele.exposure,
  eaf           = harm$eaf.exposure,          # NA in file; all 8 non-palindromic
  pval          = harm$pval.exposure,
  samplesize    = harm$samplesize.exposure,
  Phenotype     = "transferrin",
  stringsAsFactors = FALSE)
exposure_dat <- format_data(exp_raw, type="exposure")
cat("Exposure instruments:", nrow(exposure_dat), "SNPs\n")

# ---- Controls ----
controls <- list(
  list(id="ieu-a-1051", name="Transferrin saturation", role="positive control",
       expect="NEGATIVE (higher transferrin lowers saturation at fixed iron)",
       overlap="SHARES Benyamin 2014 cohort with the transferrin exposure"),
  list(id="ieu-a-89",   name="Height",                 role="negative control",
       expect="null", overlap="independent (GIANT)"),
  list(id="ieu-a-1239", name="Educational attainment", role="negative control",
       expect="null (no plausible iron pathway)", overlap="independent (SSGAC)"),
  list(id="ieu-a-1049", name="Serum iron",             role="DAG / mediator check",
       expect="tests whether iron is downstream of transferrin (mediator)",
       overlap="SHARES Benyamin 2014 cohort with the transferrin exposure")
)

rows <- list(); notes <- c()
for (ct in controls) {
  cat("\n== transferrin ->", ct$name, "(", ct$id, ") ==\n")
  out <- tryCatch(extract_outcome_data(snps=exposure_dat$SNP, outcomes=ct$id),
                  error=function(e){cat("  ERR:",conditionMessage(e),"\n");NULL})
  if (is.null(out) || nrow(out)==0) { notes<-c(notes,sprintf("- %s: no outcome data",ct$name)); next }
  dat <- harmonise_data(exposure_dat, out, action=2)
  dat <- dat[dat$mr_keep==TRUE,]
  cat("  harmonized SNPs:", nrow(dat), "\n")
  res <- mr(dat, method_list=c("mr_ivw","mr_weighted_median",
                               "mr_egger_regression","mr_weighted_mode"))
  # Steiger directionality (esp. for the iron DAG check)
  steiger <- tryCatch({
    dd <- directionality_test(dat); dd
  }, error=function(e) NULL)
  for (i in seq_len(nrow(res))) {
    rows[[length(rows)+1]] <- data.frame(
      outcome=ct$name, id=ct$id, role=ct$role, nsnp=res$nsnp[i],
      method=res$method[i], b=signif(res$b[i],4), se=signif(res$se[i],4),
      pval=signif(res$pval[i],4),
      steiger_dir=if(!is.null(steiger)) ifelse(steiger$correct_causal_direction[1],"exposure->outcome","REVERSED") else NA,
      steiger_p=if(!is.null(steiger)) signif(steiger$steiger_pval[1],3) else NA,
      overlap=ct$overlap, stringsAsFactors=FALSE)
  }
  ivw <- res[res$method=="Inverse variance weighted",]
  notes <- c(notes, sprintf("- **%s** (%s): IVW b=%.4f, p=%.3g [%d SNPs]. Expect %s. Overlap: %s.",
             ct$name, ct$role, ivw$b, ivw$pval, ivw$nsnp, ct$expect, ct$overlap))
}

res_all <- do.call(rbind, rows)
write.csv(res_all, file.path(REV,"TableS7_transferrin_controls.csv"), row.names=FALSE)

md <- c("# Table S7 — transferrin-specific controls (Part 3)","",
  sprintf("Generated %s. Exposure = 8-SNP transferrin instrument set (Benyamin ieu-a-1052), same pipeline as primary. TwoSampleMR IVW + sensitivity.",DATE_RUN),
  "TIBC positive control unavailable in OpenGWAS (Bell TIBC = deCODE GCST011368 only); near-collinear with transferrin, so omitted with note.","",
  "## Results (IVW headline per outcome)", notes, "",
  "## Interpretation",
  "- Positive control (transferrin saturation): a strong NEGATIVE estimate validates the transferrin instruments' biology (more transferrin -> more unsaturated binding sites -> lower saturation), i.e. the NTBI-buffering premise. NOTE: TSAT shares the Benyamin cohort, so this is a within-sample biological-consistency check, not out-of-sample.",
  "- Negative controls (height, educational attainment): null estimates confirm the transferrin instruments are not pleiotropically associated with unrelated traits.",
  "- DAG / mediator check (transferrin -> serum iron): a null or weak estimate supports treating serum iron as a PARALLEL exposure (conditioned on in MVMR), not a mediator on the transferrin->tau path; a strong estimate would raise a mediator concern. Steiger direction reported. Report the result honestly either way.")
writeLines(md, file.path(REV,"TableS7_transferrin_controls_summary.md"))

cat("\n=== TableS7 ===\n"); print(res_all[,c("outcome","role","nsnp","method","b","se","pval","steiger_dir")], row.names=FALSE)
