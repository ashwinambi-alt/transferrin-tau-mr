##############################################################################
# R1_G_phenoscanner.R
#
# Reviewer response G: PhenoScanner v2 pleiotropy lookup for the 8 primary
# transferrin instruments. Per-SNP loop with retries because the batch
# endpoint is flaky.
#
# Writes: results/reviewer_response/G_phenoscanner_raw.csv
#         results/reviewer_response/G_phenoscanner_summary.csv
##############################################################################

args <- commandArgs(trailingOnly = TRUE)
PROJECT_DIR <- if (length(args) >= 1) args[1] else normalizePath("..", mustWork = TRUE)
OUT_DIR     <- file.path(PROJECT_DIR, "results", "reviewer_response")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(phenoscanner)
  library(dplyr)
})

options(timeout = 240)

snps <- c("rs1495741", "rs174577", "rs1800562", "rs744653",
          "rs7646473", "rs8177240", "rs9268633", "rs9990333")

query_one <- function(rsid, tries = 3, sleep = 6) {
  for (k in seq_len(tries)) {
    cat(sprintf("  %s (attempt %d)...\n", rsid, k))
    r <- tryCatch(
      phenoscanner(snpquery = rsid, catalogue = "GWAS",
                   pvalue = 1e-5, proxies = "EUR",
                   r2 = 0.8, build = 37),
      error = function(e) { cat("    ERR:", conditionMessage(e), "\n"); NULL }
    )
    if (!is.null(r) && !is.null(r$results)) {
      cat(sprintf("    OK: %d rows\n", nrow(r$results)))
      return(r)
    }
    Sys.sleep(sleep)
  }
  cat("    GIVING UP on", rsid, "\n")
  NULL
}

all_results <- list(); all_meta <- list()
for (rsid in snps) {
  r <- query_one(rsid)
  if (!is.null(r)) {
    if (nrow(r$results) > 0) all_results[[rsid]] <- r$results
    if (!is.null(r$snps) && nrow(r$snps) > 0) all_meta[[rsid]] <- r$snps
  }
  Sys.sleep(3)   # be polite between SNPs
}

if (length(all_results) == 0) {
  cat("\nPhenoScanner returned nothing — fallback to GWAS Catalog REST.\n")
  library(httr); library(jsonlite)
  gc_query <- function(rsid) {
    url <- sprintf("https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/%s/associations", rsid)
    r <- tryCatch(GET(url, timeout(60)), error = function(e) NULL)
    if (is.null(r) || http_error(r)) return(NULL)
    j <- content(r, as = "parsed", type = "application/json")
    if (is.null(j$`_embedded`$associations)) return(NULL)
    do.call(rbind, lapply(j$`_embedded`$associations, function(a) {
      data.frame(
        snp   = rsid,
        pvalue = ifelse(is.null(a$pvalue), NA, a$pvalue),
        beta   = ifelse(is.null(a$betaNum), NA, a$betaNum),
        trait  = paste(sapply(a$`_links`$efoTraits$href, function(x) x), collapse=";"),
        stringsAsFactors = FALSE
      )
    }))
  }
  gc_all <- do.call(rbind, lapply(snps, function(rs) {
    cat("  GWAS Catalog", rs, "...\n"); Sys.sleep(2); gc_query(rs)
  }))
  if (!is.null(gc_all) && nrow(gc_all) > 0) {
    write.csv(gc_all, file.path(OUT_DIR, "G_gwas_catalog_fallback.csv"), row.names = FALSE)
    cat("Wrote GWAS Catalog fallback.\n")
  }
  quit(status = 0)
}

results_df <- do.call(rbind, all_results)
snpmeta_df <- if (length(all_meta) > 0) do.call(rbind, all_meta) else NULL

write.csv(results_df, file.path(OUT_DIR, "G_phenoscanner_raw.csv"), row.names = FALSE)
if (!is.null(snpmeta_df)) {
  write.csv(snpmeta_df, file.path(OUT_DIR, "G_phenoscanner_snpmeta.csv"), row.names = FALSE)
}
cat(sprintf("\nTotal rows: %d across %d SNPs with hits\n",
            nrow(results_df), length(unique(results_df$snp))))

# ---- categorise ----
iron_terms     <- "(?i)iron|ferritin|transferrin|haemoglobin|hemoglobin|MCH|MCV|red cell|erythrocyte|hepcidin|haematocrit|hematocrit|reticulocyte"
tau_terms      <- "(?i)tau|neurofibrillary|alzheimer|dementia|cognitive"
amyloid_terms  <- "(?i)amyloid|abeta|a-beta"
lipid_terms    <- "(?i)cholesterol|LDL|HDL|triglyc|lipid|fatty acid|phospholipid|sphingo|apolipo"
immune_terms   <- "(?i)immune|autoimmune|inflammat|rheumatoid|lupus|celiac|coeliac|IgG|IgA|IgE|hepatitis|HLA|MHC|psoriasis|multiple sclerosis|IBD|Crohn|ulcerative"
metabolic_terms<- "(?i)diabetes|glucose|insulin|BMI|adiposity|obes|body mass|waist|hip"
liver_terms    <- "(?i)liver|hepatic|ALT|AST|GGT|bilirubin"

results_df <- results_df %>%
  mutate(
    cat_iron    = grepl(iron_terms, trait),
    cat_tau     = grepl(tau_terms, trait),
    cat_amyloid = grepl(amyloid_terms, trait),
    cat_lipid   = grepl(lipid_terms, trait),
    cat_immune  = grepl(immune_terms, trait),
    cat_metab   = grepl(metabolic_terms, trait),
    cat_liver   = grepl(liver_terms, trait)
  )

summary_df <- results_df %>%
  group_by(snp) %>%
  summarise(
    n_hits      = n(),
    n_iron      = sum(cat_iron),
    n_tau_dem   = sum(cat_tau),
    n_amyloid   = sum(cat_amyloid),
    n_lipid     = sum(cat_lipid),
    n_immune    = sum(cat_immune),
    n_metabolic = sum(cat_metab),
    n_liver     = sum(cat_liver),
    n_other     = sum(!cat_iron & !cat_tau & !cat_amyloid & !cat_lipid &
                        !cat_immune & !cat_metab & !cat_liver),
    top_non_iron_traits =
      paste(head(unique(trait[!cat_iron]), 8), collapse = "; "),
    .groups = "drop"
  )

full <- data.frame(snp = snps, stringsAsFactors = FALSE)
summary_df <- merge(full, summary_df, by = "snp", all.x = TRUE)
num_cols <- sapply(summary_df, is.numeric)
summary_df[, num_cols][is.na(summary_df[, num_cols])] <- 0

write.csv(summary_df, file.path(OUT_DIR, "G_phenoscanner_summary.csv"),
          row.names = FALSE)

cat("\n== Analysis G: PhenoScanner summary ==\n")
print(summary_df[, c("snp", "n_hits", "n_iron", "n_tau_dem", "n_amyloid",
                     "n_lipid", "n_immune", "n_metabolic", "n_liver", "n_other")])
