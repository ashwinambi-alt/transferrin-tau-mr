##############################################################################
# R1_G_gwas_catalog.R
#
# Reviewer response G: pleiotropy lookup for all 8 primary transferrin
# instruments using the GWAS Catalog REST API (PhenoScanner is currently
# unreachable).
#
# For each SNP:
#   1. GET  /singleNucleotidePolymorphisms/{rsid}/associations
#   2. For each association, follow _links.study.href to get the
#      disease/trait label and p-value threshold.
# Then categorise each trait into iron / tau / amyloid / lipid / immune /
# metabolic / liver / other, matching the categories the reviewer
# specifically flagged.
##############################################################################

args <- commandArgs(trailingOnly = TRUE)
PROJECT_DIR <- if (length(args) >= 1) args[1] else normalizePath("..", mustWork = TRUE)
OUT_DIR     <- file.path(PROJECT_DIR, "results", "reviewer_response")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(httr); library(jsonlite); library(dplyr)
})

options(timeout = 120)

snps <- c("rs1495741", "rs174577", "rs1800562", "rs744653",
          "rs7646473", "rs8177240", "rs9268633", "rs9990333")

# Manual gene annotation (from manuscript Supp Table S1 and dbSNP)
gene_map <- c(
  rs1495741 = "NAT2",   rs174577  = "FADS1/FADS2", rs1800562 = "HFE",
  rs744653  = "TFRC",   rs7646473 = "TF (5' region)",
  rs8177240 = "TF",     rs9268633 = "HLA-DRB1/DQA1", rs9990333 = "TFRC downstream"
)

get_json <- function(url) {
  r <- tryCatch(GET(url, timeout(60)), error = function(e) NULL)
  if (is.null(r) || http_error(r)) return(NULL)
  content(r, as = "parsed", type = "application/json", simplifyVector = FALSE)
}

fetch_assocs <- function(rsid) {
  cat(sprintf("  %s ...", rsid))
  url <- sprintf("https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/%s/associations", rsid)
  root <- get_json(url)
  if (is.null(root) || is.null(root$`_embedded`$associations)) {
    cat(" no data\n"); return(NULL)
  }
  assocs <- root$`_embedded`$associations
  cat(sprintf(" %d assocs\n", length(assocs)))
  out <- do.call(rbind, lapply(assocs, function(a) {
    study_url <- a$`_links`$study$href
    study <- if (!is.null(study_url)) get_json(study_url) else NULL
    trait <- if (!is.null(study$diseaseTrait$trait)) study$diseaseTrait$trait else NA
    pubmed <- if (!is.null(study$publicationInfo$pubmedId)) study$publicationInfo$pubmedId else NA
    n_case  <- tryCatch(sum(unlist(lapply(study$ancestries,
                function(x) if (!is.null(x$numberOfIndividuals)) x$numberOfIndividuals else 0))), error=function(e) NA)
    pval <- if (!is.null(a$pvalue)) a$pvalue else NA
    or   <- if (!is.null(a$orPerCopyNum)) a$orPerCopyNum else NA
    beta <- if (!is.null(a$betaNum)) a$betaNum else NA
    beta_unit <- if (!is.null(a$betaUnit)) a$betaUnit else NA
    data.frame(
      snp    = rsid,
      trait  = trait,
      pvalue = pval,
      beta   = beta,
      beta_unit = beta_unit,
      or     = or,
      pubmed = pubmed,
      n_reported = n_case,
      stringsAsFactors = FALSE
    )
  }))
  Sys.sleep(1)
  out
}

all <- do.call(rbind, lapply(snps, fetch_assocs))
if (is.null(all)) stop("No GWAS Catalog data returned for any SNP")

all$gene <- gene_map[all$snp]

write.csv(all, file.path(OUT_DIR, "G_gwas_catalog_raw.csv"), row.names = FALSE)
cat(sprintf("\nTotal associations: %d across %d SNPs\n",
            nrow(all), length(unique(all$snp))))

# ---- categorise traits ----
iron_terms     <- "(?i)iron|ferritin|transferrin|haemoglobin|hemoglobin|MCH\\b|MCV\\b|mean corpuscular|red cell|red blood|erythrocyte|hepcidin|haematocrit|hematocrit|reticulocyte"
tau_terms      <- "(?i)\\btau\\b|neurofibrillary|alzheimer|dementia|cognitive"
amyloid_terms  <- "(?i)amyloid|abeta|a-beta|APP\\b"
lipid_terms    <- "(?i)cholesterol|\\bLDL\\b|\\bHDL\\b|triglyc|lipid|fatty acid|phospholipid|sphingo|apolipo"
immune_terms   <- "(?i)immune|autoimmune|inflammat|rheumatoid|lupus|celiac|coeliac|IgG|IgA|IgE|hepatitis|psoriasis|multiple sclerosis|inflammatory bowel|Crohn|ulcerative|type 1 diabetes|asthma|allerg"
metabolic_terms<- "(?i)type 2 diabetes|glucose|insulin|\\bBMI\\b|adiposity|obes|body mass|waist|hip circumference"
liver_terms    <- "(?i)liver|hepatic|\\bALT\\b|\\bAST\\b|\\bGGT\\b|bilirubin|non-alcoholic fatty|NAFLD"

all$trait <- ifelse(is.na(all$trait), "", all$trait)
all <- all %>%
  mutate(
    cat_iron    = grepl(iron_terms, trait),
    cat_tau     = grepl(tau_terms, trait),
    cat_amyloid = grepl(amyloid_terms, trait),
    cat_lipid   = grepl(lipid_terms, trait),
    cat_immune  = grepl(immune_terms, trait),
    cat_metab   = grepl(metabolic_terms, trait),
    cat_liver   = grepl(liver_terms, trait)
  )

summary_df <- all %>%
  group_by(snp) %>%
  summarise(
    gene = gene_map[first(snp)],
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
    non_iron_traits =
      paste(head(unique(trait[!cat_iron & trait != ""]), 12), collapse = "; "),
    .groups = "drop"
  )
full <- data.frame(snp = snps, gene = gene_map[snps], stringsAsFactors = FALSE)
summary_df <- merge(full, summary_df, by = c("snp","gene"), all.x = TRUE)
num_cols <- sapply(summary_df, is.numeric)
summary_df[, num_cols][is.na(summary_df[, num_cols])] <- 0

# concern flag
summary_df$pleiotropy_concern <- with(summary_df,
  ifelse(n_amyloid + n_tau_dem > 0, "HIGH — direct outcome-adjacent trait",
  ifelse(n_lipid + n_immune + n_liver >= 2, "MEDIUM — multiple non-iron systems",
  ifelse(n_other > n_iron & n_other >= 3,   "MEDIUM — many non-iron hits",
  ifelse(n_other + n_lipid + n_immune + n_liver + n_metabolic == 0, "LOW — iron-only",
                                            "LOW — few non-iron hits")))))

write.csv(summary_df, file.path(OUT_DIR, "G_gwas_catalog_summary.csv"),
          row.names = FALSE)

cat("\n== Analysis G: GWAS Catalog pleiotropy summary ==\n")
print(summary_df[, c("snp","gene","n_hits","n_iron","n_tau_dem","n_amyloid",
                     "n_lipid","n_immune","n_metabolic","n_liver","n_other",
                     "pleiotropy_concern")])
