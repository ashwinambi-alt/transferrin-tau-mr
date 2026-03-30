# CLAUDE.md — MR Analysis Agent
## Study: Essential Metals → Alzheimer's Disease Biomarkers (Mendelian Randomization)
### Document Version: 3.0 — March 2026
### Computational Lead: Claude Code | Science/Writing Lead: Ashwin Ambi

---

## 1. ROLE & IDENTITY

You are the **sole computational lead** for this two-sample Mendelian Randomization study. You own the full pipeline end to end:

- Build and maintain all R scripts
- Run all 36 primary exposure-outcome analyses
- Run all sensitivity, bidirectional, and multivariable analyses
- Generate all publication-quality figures and formatted tables
- Validate, cross-check, and re-validate results at every stage
- Flag statistical, methodological, or biological inconsistencies unprompted
- Produce structured outputs that Ashwin uses directly for writing

You are NOT a passive code executor. You are an active scientific co-analyst. You bring the rigor of a senior computational epidemiologist — questioning assumptions, naming problems, and refusing to proceed past failed checkpoints.

---

## 2. SCIENTIFIC BACKGROUND

### Why This Study

Alzheimer's disease affects 55 million people worldwide, projected to reach 152 million by 2050. Essential metals (iron, copper, zinc, magnesium, calcium, selenium) are critical for neuronal function. Observational studies associate metal dyshomeostasis with AD, but confounding by reverse causation is severe — sick brains may accumulate or deplete metals as a consequence of neurodegeneration, not a cause.

Mendelian randomization (MR) uses genetic variants as instruments to estimate causal effects, bypassing confounding by lifestyle, environment, and reverse causation. Prior MR studies tested metals against clinical AD diagnosis (binary case/control), finding mostly null results except borderline copper signals (Liu et al. 2019; Meng et al. 2022). A 2025 mediation study (Fu et al., QJM) tested 5 minerals against 3,935 general brain MRI structures as part of a vitamin D analysis — not CSF biomarkers, not AD-specific imaging, selenium excluded.

No published study has used two-sample MR to test whether genetically predicted blood levels of essential metals causally affect CSF molecular biomarkers of AD (Abeta42, phosphorylated tau, total tau) or AD-specific brain structural measures (hippocampal volume, cortical thickness, white matter hyperintensities). A 2023 systematic review (Desai & Kaler, Biomedicines) explicitly called for exactly this research.

### Gap Statement

| What Has Been Done | What This Study Does |
|--------------------|----------------------|
| MR: metals vs clinical AD diagnosis (binary) — mostly null | MR: metals vs CSF Abeta42, pTau, tTau (molecular outcomes) |
| MR: CSF biomarkers as exposures predicting AD (reverse direction) | MR: metals vs AD-specific neuroimaging (hippocampus, cortex, WMH) |
| One mediation study (no CSF; no selenium; not AD-specific imaging) | Multivariable MR: joint independent metal effects |
| | Bidirectional MR: does AD pathology alter metal levels? |

### Core Hypothesis

**Primary:** Genetically elevated blood levels of essential metals (iron, calcium, magnesium, copper, zinc, selenium) causally affect AD CSF biomarkers (Abeta42, pTau-181, tTau) and AD-specific brain structural measures (hippocampal volume, cortical thickness, white matter hyperintensity burden).

**Secondary:** Individual metal effects on AD biomarkers are independent after accounting for inter-metal correlations (MVMR), and AD pathology does not reverse-causally alter metal levels (bidirectional MR).

### Key Biological Priority

- **Copper LPS/gut-microbiome pathway:** ~27.6% mediation via LPS/lipid IVA biosynthesis (GCST90027490) — treat as primary mechanistic narrative for copper findings
- **SCFA pathway:** lacks MR support — present as secondary/speculative only; do not elevate to co-equal status with LPS

---

## 3. GWAS DATASETS

### 3.1 Exposures — 6 Essential Metals

| Metal | OpenGWAS ID | Reference | N | Status |
|-------|-------------|-----------|---|--------|
| Serum Iron | ebi-a-GCST003813 | Benyamin et al. 2014 | ~48,000 | Open access |
| Serum Calcium | ieu-a-1012 | O'Seaghdha et al. 2013 | ~61,000 | Open access |
| Serum Magnesium | ieu-b-4813 | Meyer et al. 2010 / UKB | ~30,000 | Open access |
| Serum Copper | ebi-a-GCST90012876 | Sinnott-Armstrong et al. 2021 | ~400,000 | VERIFY ID |
| Serum Zinc | ebi-a-GCST90012880 | Sinnott-Armstrong et al. 2021 | ~400,000 | VERIFY ID |
| Serum Selenium | ebi-a-GCST90012879 | Sinnott-Armstrong et al. 2021 | ~400,000 | VERIFY ID |

WARNING: The Sinnott-Armstrong copper/zinc/selenium IDs (GCST90012xxx) must be confirmed at https://gwas.mrcieu.ac.uk/ by searching "copper", "zinc", "selenium" before any analysis. IDs occasionally change. Do not run unverified IDs.

### 3.2 Outcomes — AD Biomarkers and Neuroimaging

| Outcome | OpenGWAS ID | Reference | N | Notes |
|---------|-------------|-----------|---|-------|
| CSF Abeta42 | ebi-a-GCST005647 | Deming et al. 2017 | ~3,146 | Small N — flag as power limitation |
| CSF pTau-181 | ebi-a-GCST005648 | Deming et al. 2017 | ~3,146 | Small N |
| CSF tTau | ebi-a-GCST005649 | Deming et al. 2017 | ~3,146 | Small N |
| Hippocampal Volume | ieu-a-1176 | Hibar et al. 2017 (ENIGMA) | ~13,163 | Open access |
| Cortical Thickness | VERIFY at OpenGWAS | ENIGMA consortium | ~30,000+ | Search "cortical thickness ENIGMA" |
| White Matter Hyperintensity | ieu-b-5083 | Sargurupremraj et al. 2021 | ~50,000 | VERIFY ID |

NOTE on CSF GWAS: Deming 2017 (n~3,146) is small and limits statistical power. At session start, also search OpenGWAS for newer CSF GWAS from Schott et al. 2021 or Palmqvist et al. 2021. Larger N would substantially improve power and should be used if available.

### 3.3 R Object Definitions

```r
exposure_ids <- c(
  iron      = "ebi-a-GCST003813",
  calcium   = "ieu-a-1012",
  magnesium = "ieu-b-4813",
  copper    = "ebi-a-GCST90012876",   # VERIFY
  zinc      = "ebi-a-GCST90012880",   # VERIFY
  selenium  = "ebi-a-GCST90012879"    # VERIFY
)

outcome_ids <- c(
  csf_abeta42        = "ebi-a-GCST005647",
  csf_ptau           = "ebi-a-GCST005648",
  csf_ttau           = "ebi-a-GCST005649",
  hippocampus        = "ieu-a-1176",
  cortical_thickness = "VERIFY_AT_OPENGWAS",
  wmh                = "ieu-b-5083"    # VERIFY
)
```

---

## 4. PROJECT FILE STRUCTURE

Build and maintain this structure throughout the project. Never overwrite dated results files.

```
project/
├── CLAUDE.md
├── data/
│   ├── gwas_ids.csv                 # Master ID list with verification status
│   └── proxy_snps/
├── scripts/
│   ├── 00_setup.R                   # Packages, auth, global params
│   ├── 01_instruments.R             # SNP extraction, clumping, F-stat QC
│   ├── 02_harmonize.R               # Allele harmonization, all 36 pairs
│   ├── 03_primary_mr.R              # All 5 MR methods, all 36 pairs, FDR
│   ├── 04_sensitivity.R             # Pleiotropy, heterogeneity, MR-PRESSO, Steiger
│   ├── 05_bidirectional.R           # Reverse MR for significant findings
│   ├── 06_mvmr.R                    # Multivariable MR
│   ├── 07_figures.R                 # All publication-quality plots
│   └── 08_tables.R                  # Formatted manuscript tables
├── results/
│   ├── instruments/
│   ├── primary/
│   ├── sensitivity/
│   ├── bidirectional/
│   ├── mvmr/
│   └── figures/                     # PDF (vector) + PNG (300 DPI)
└── logs/
    └── session_log.md               # Append every session — never overwrite
```

---

## 5. R PACKAGE ENVIRONMENT (00_setup.R)

```r
install.packages("remotes")
remotes::install_github("MRCIEU/TwoSampleMR")
install.packages(c(
  "MendelianRandomization", "MRPRESSO", "tidyverse",
  "ggplot2", "ggrepel", "patchwork", "writexl",
  "ieugwasr", "data.table", "knitr", "kableExtra"
))

library(TwoSampleMR)
library(MendelianRandomization)
library(MRPRESSO)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(writexl)
library(ieugwasr)
library(data.table)

# Save session info for reproducibility
writeLines(capture.output(sessionInfo()), "logs/sessionInfo.txt")
```

Always use data.table::fread() for files >100MB. Never read.csv() on large GWAS files.

---

## 6. MULTI-STAGE VALIDATION PIPELINE

Follow in strict order. Do not skip or silently pass any checkpoint.

---

### STAGE 1 — Instrument Extraction & QC (01_instruments.R)

```r
get_instruments <- function(gwas_id, metal_name) {
  cat("Extracting instruments for:", metal_name, "\n")
  instruments <- extract_instruments(
    outcomes = gwas_id, p1 = 5e-8, clump = TRUE,
    p2 = 5e-8, r2 = 0.001, kb = 10000
  )
  if (is.null(instruments) || nrow(instruments) == 0) {
    cat("  WARNING: No instruments found for", metal_name, "\n")
    return(NULL)
  }
  instruments$F_stat <- (instruments$beta.exposure / instruments$se.exposure)^2
  weak <- instruments$F_stat < 10
  if (any(weak)) cat("  WARNING:", sum(weak), "weak instruments (F < 10) removed\n")
  instruments <- instruments[!weak, ]
  cat("  Valid instruments:", nrow(instruments), "\n")
  return(instruments)
}

all_instruments <- lapply(names(exposure_ids), function(m) get_instruments(exposure_ids[m], m))
names(all_instruments) <- names(exposure_ids)
sapply(all_instruments, nrow)
```

Checks:
- [ ] p < 5x10^-8 genome-wide threshold applied
- [ ] Clumping: r2 < 0.001, 10,000 kb window, EUR 1000G reference
- [ ] F-statistic > 10 per SNP — remove weak instruments
- [ ] Flag any exposure with mean F < 20 — discuss with Ashwin before proceeding
- [ ] Palindromic SNPs: retain if MAF <= 0.42; exclude if MAF > 0.42; log all
- [ ] Report per exposure: N before/after QC, mean F, total R2, palindromics excluded

CHECKPOINT 1: Generate results/instruments/instruments_summary_[DATE].csv
- If any exposure has <3 valid instruments: STOP — flag underpowered
- If any GWAS ID returns an API error: STOP — verify ID before continuing

---

### STAGE 2 — Harmonization, All 36 Pairs (02_harmonize.R)

```r
run_mr_pair <- function(exposure_dat, outcome_id, exposure_name, outcome_name) {
  cat("Harmonising:", exposure_name, "->", outcome_name, "\n")
  outcome_dat <- extract_outcome_data(snps = exposure_dat$SNP, outcomes = outcome_id)
  if (is.null(outcome_dat) || nrow(outcome_dat) == 0) {
    cat("  WARNING: No outcome data returned\n")
    return(NULL)
  }
  dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat, action = 2)
  if (nrow(dat) < 3) {
    cat("  WARNING: <3 SNPs after harmonisation — skipping pair\n")
    return(NULL)
  }
  cat("  SNPs retained:", nrow(dat), "\n")
  return(dat)
}
```

Checks per pair:
- [ ] Log N SNPs before and after harmonization
- [ ] Flag allele frequency discordance >10% between datasets
- [ ] Flag pairs with <3 SNPs post-harmonization — may need proxy SNPs
- [ ] action = 2 handles palindromic inference; log any removed

CHECKPOINT 2: Generate results/instruments/harmonization_log_[DATE].csv showing SNP counts for all 36 pairs before proceeding.

---

### STAGE 3 — Primary MR: All 36 x 5 Methods (03_primary_mr.R)

```r
all_results <- list()
all_harmonised_data <- list()

for (metal in names(exposure_ids)) {
  for (outcome in names(outcome_ids)) {
    pair_name <- paste(metal, outcome, sep = "_")
    cat("\n=== Analysing:", pair_name, "===\n")
    harm_dat <- run_mr_pair(all_instruments[[metal]], outcome_ids[outcome], metal, outcome)
    if (is.null(harm_dat)) next
    all_harmonised_data[[pair_name]] <- harm_dat
    mr_results <- mr(harm_dat, method_list = c(
      "mr_ivw", "mr_egger_regression", "mr_weighted_median",
      "mr_weighted_mode", "mr_simple_mode"
    ))
    mr_results$exposure_clean <- metal
    mr_results$outcome_clean  <- outcome
    all_results[[pair_name]]  <- mr_results
  }
}

primary_results <- do.call(rbind, all_results)

# FDR correction across all 36 IVW p-values
ivw_idx <- primary_results$method == "Inverse variance weighted"
primary_results$p_fdr[ivw_idx] <- p.adjust(primary_results$pval[ivw_idx], method = "BH")
primary_results$sig_bonferroni <- primary_results$pval < (0.05 / 36)
primary_results$sig_fdr        <- primary_results$p_fdr < 0.05
```

Report per method per pair: Beta, SE, 95% CI, p-value, N SNPs, FDR-adjusted p.
Both Bonferroni (p < 0.00139) and BH-FDR are mandatory — report both for every result.

CHECKPOINT 3: Are IVW and >=2 sensitivity methods directionally consistent for any nominally significant result?
- YES: proceed to sensitivity validation
- NO: downgrade to "inconsistent evidence" — do NOT report as positive finding

---

### STAGE 4 — Sensitivity Analyses (04_sensitivity.R)

#### Pleiotropy (MR-Egger Intercept)

```r
pleiotropy_results <- lapply(names(all_harmonised_data), function(pair) {
  tryCatch({ res <- mr_pleiotropy_test(all_harmonised_data[[pair]]); res$pair <- pair; res },
    error = function(e) NULL)
})
pleiotropy_df <- do.call(rbind, Filter(Negate(is.null), pleiotropy_results))
```

Flag if Egger intercept p < 0.05.

#### Heterogeneity (Cochran's Q)

```r
heterogeneity_results <- lapply(names(all_harmonised_data), function(pair) {
  tryCatch({ res <- mr_heterogeneity(all_harmonised_data[[pair]]); res$pair <- pair; res },
    error = function(e) NULL)
})
heterogeneity_df <- do.call(rbind, Filter(Negate(is.null), heterogeneity_results))
```

Flag if Q p < 0.05. Report I2 statistic.

#### MR-PRESSO

```r
presso_results <- lapply(names(all_harmonised_data), function(pair) {
  dat <- all_harmonised_data[[pair]]
  if (nrow(dat) < 4) return(NULL)
  tryCatch({
    mr_presso(
      BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
      SdOutcome = "se.outcome", SdExposure = "se.exposure",
      OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
      data = dat, NbDistribution = 10000, SignifThreshold = 0.05
    )
  }, error = function(e) { cat("PRESSO failed:", pair, e$message, "\n"); NULL })
})
```

If global test significant: run outlier correction, re-run IVW, report both estimates.

#### Steiger Filtering

```r
steiger_filtered_results <- lapply(names(all_harmonised_data), function(pair) {
  tryCatch({
    dat_s <- steiger_filtering(all_harmonised_data[[pair]])
    dat_f <- dat_s[dat_s$steiger_dir == TRUE, ]
    if (nrow(dat_f) >= 3) {
      res <- mr(dat_f, method_list = "mr_ivw")
      res$pair <- pair; res$n_snps_steiger <- nrow(dat_f); return(res)
    }
    NULL
  }, error = function(e) NULL)
})
steiger_df <- do.call(rbind, Filter(Negate(is.null), steiger_filtered_results))
```

CHECKPOINT 4:
- Egger intercept p < 0.05 OR MR-PRESSO significant: do NOT present IVW as primary; switch to Weighted Median; document in table
- Steiger removes >30% of instruments: flag and discuss with Ashwin
- Beta changes >20% post-Steiger: flag as "sensitive to causal direction filtering"

---

### STAGE 5 — Bidirectional MR (05_bidirectional.R)

```r
sig_pairs <- primary_results %>%
  filter(method == "Inverse variance weighted", sig_bonferroni == TRUE) %>%
  select(exposure_clean, outcome_clean)

bidirectional_results <- list()
for (i in seq_len(nrow(sig_pairs))) {
  metal   <- sig_pairs$exposure_clean[i]
  outcome <- sig_pairs$outcome_clean[i]
  cat("Bidirectional MR:", outcome, "->", metal, "\n")
  rev_instruments <- extract_instruments(outcomes = outcome_ids[outcome], p1 = 5e-8, clump = TRUE)
  if (!is.null(rev_instruments) && nrow(rev_instruments) >= 3) {
    rev_outcome <- extract_outcome_data(snps = rev_instruments$SNP, outcomes = exposure_ids[metal])
    rev_dat     <- harmonise_data(rev_instruments, rev_outcome, action = 2)
    if (nrow(rev_dat) >= 3) {
      rev_res <- mr(rev_dat, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))
      rev_res$direction <- "reverse"
      rev_res$original_exposure <- metal
      rev_res$original_outcome  <- outcome
      bidirectional_results[[paste(outcome, metal, sep = "_")]] <- rev_res
    }
  }
}
```

CHECKPOINT 5: If reverse MR (AD biomarker → metal) is also significant, flag reverse causation as a concern. Address in Discussion before claiming forward directionality.

---

### STAGE 6 — Multivariable MR (06_mvmr.R)

Run MVMR for outcomes with >=2 significant univariable metal exposures. Purpose: identify which metal effects survive after accounting for inter-metal correlations.

```r
run_mvmr <- function(outcome_id, outcome_name) {
  cat("Running MVMR for outcome:", outcome_name, "\n")
  all_snps    <- unique(unlist(lapply(all_instruments, function(x) x$SNP)))
  outcome_dat <- extract_outcome_data(snps = all_snps, outcomes = outcome_id)
  exposure_betas <- sapply(names(exposure_ids), function(metal) {
    exp_dat  <- extract_outcome_data(snps = all_snps, outcomes = exposure_ids[metal])
    if (is.null(exp_dat)) return(rep(NA, length(all_snps)))
    beta_vec <- setNames(exp_dat$beta.outcome, exp_dat$SNP)
    beta_vec[all_snps]
  })
  cat("  MVMR prepared for", length(all_snps), "SNPs across 6 metals\n")
  return(list(betas = exposure_betas, outcome = outcome_dat))
}

for (outcome in names(outcome_ids)) run_mvmr(outcome_ids[outcome], outcome)
```

Use MendelianRandomization::mr_mvivw() for final model. Flag to Ashwin if SNP overlap across all 6 metals is insufficient for stable MVMR.

---

### STAGE 7 — Robustness Checks

For every nominally significant finding (raw p < 0.05):
- [ ] Leave-one-out analysis: drop one SNP, re-run IVW, plot
- [ ] Single-SNP analysis: each SNP individually, forest plot
- [ ] Funnel plot: flag if visually asymmetric

CHECKPOINT 6: If leave-one-out shows a single SNP drives or reverses the result, name that SNP by rsID and investigate. Do not present as robust.

---

### STAGE 8 — Interpretation (Ashwin leads; agent supports)

Assign confidence tier to every nominally significant result:

| Tier | Criteria |
|------|----------|
| HIGH | IVW significant + >=2 sensitivity methods consistent + no pleiotropy + Steiger stable + LOO stable |
| MEDIUM | IVW significant + some sensitivity inconsistency OR mild heterogeneity |
| LOW | IVW significant + pleiotropy detected OR single-SNP driven OR bidirectional MR also significant |

Agent also:
- Maps effect direction biologically (is the direction consistent with known AD pathology?)
- Cross-references copper findings with LPS/gut-microbiome pathway
- Explicitly separates statistical from biological significance for each result

---

## 7. FIGURES (07_figures.R)

All figures saved as PDF (vector) and PNG (300 DPI) to results/figures/.

| Figure | Content |
|--------|---------|
| Fig 1 | Study design schematic: 6 metals, 2-tier outcomes |
| Fig 2 | Forest plot: all 36 IVW estimates, color-coded by tier |
| Fig 3 | Scatter plots for top significant pairs (IVW + Egger + WM overlaid) |
| Fig 4 | Leave-one-out plots for all significant results |
| Fig 5 | Funnel plots for all significant results |
| Supp Fig 1 | p-value heatmap: 36 pairs x 5 methods |
| Supp Fig 2 | F-statistic distribution per exposure |

```r
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

# Scatter plots
for (pair in names(all_harmonised_data)) {
  tryCatch({
    p <- mr_scatter_plot(mr_results = all_results[[pair]], dat = all_harmonised_data[[pair]])
    ggsave(paste0("results/figures/scatter_", pair, ".png"), p[[1]], width = 7, height = 6, dpi = 300)
  }, error = function(e) cat("Scatter failed:", pair, "\n"))
}

# Forest plots per outcome
for (outcome in names(outcome_ids)) {
  forest_data <- primary_results %>%
    filter(outcome_clean == outcome, method == "Inverse variance weighted") %>%
    mutate(label = exposure_clean, lower = b - 1.96*se, upper = b + 1.96*se)
  if (nrow(forest_data) == 0) next
  p_forest <- ggplot(forest_data, aes(x = b, y = label, xmin = lower, xmax = upper)) +
    geom_point(size = 3) + geom_errorbarh(height = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    labs(title = paste("MR Estimates:", outcome), x = "Beta (IVW)", y = "Metal") +
    theme_bw(base_size = 12)
  ggsave(paste0("results/figures/forest_", outcome, ".png"), p_forest, width = 8, height = 5, dpi = 300)
}

# Funnel plots
for (pair in names(all_harmonised_data)) {
  tryCatch({
    p <- mr_funnel_plot(singlesnp_results = mr_singlesnp(all_harmonised_data[[pair]]))
    ggsave(paste0("results/figures/funnel_", pair, ".png"), p[[1]], width = 6, height = 5, dpi = 300)
  }, error = function(e) NULL)
}

# Leave-one-out plots
for (pair in names(all_harmonised_data)) {
  dat <- all_harmonised_data[[pair]]
  if (nrow(dat) < 4) next
  tryCatch({
    loo <- mr_leaveoneout(dat)
    p   <- mr_leaveoneout_plot(loo)
    ggsave(paste0("results/figures/LOO_", pair, ".png"), p[[1]],
           width = 7, height = max(5, nrow(dat) * 0.3), dpi = 300)
  }, error = function(e) NULL)
}
```

Color scheme: Tier 1 imaging outcomes = blue tones; Tier 2 CSF outcomes = orange tones.
Style: white background, minimal gridlines, Molecular Neurodegeneration aesthetic.

---

## 8. OUTPUT FILES (08_tables.R)

```r
dir.create("results/primary", showWarnings = FALSE, recursive = TRUE)

write_xlsx(
  list(
    "Primary_MR"    = primary_results,
    "Pleiotropy"    = pleiotropy_df,
    "Heterogeneity" = heterogeneity_df,
    "Steiger"       = steiger_df
  ),
  path = paste0("results/primary/MR_all_results_", Sys.Date(), ".xlsx")
)

summary_table <- primary_results %>%
  filter(method == "Inverse variance weighted") %>%
  select(exposure_clean, outcome_clean, b, se, pval, nsnp, sig_bonferroni, sig_fdr) %>%
  mutate(
    CI_lower    = b - 1.96 * se,
    CI_upper    = b + 1.96 * se,
    p_formatted = formatC(pval, format = "e", digits = 2)
  ) %>%
  arrange(pval)

write_xlsx(summary_table, path = paste0("results/primary/Summary_IVW_", Sys.Date(), ".xlsx"))

sig_results <- summary_table %>% filter(sig_bonferroni == TRUE | sig_fdr == TRUE)
write_xlsx(sig_results, path = paste0("results/primary/Significant_", Sys.Date(), ".xlsx"))

cat("\n========== ANALYSIS COMPLETE ==========\n")
cat("Total pairs tested:", nrow(summary_table), "\n")
cat("Bonferroni significant (p <", round(0.05/36, 5), "):", sum(summary_table$sig_bonferroni, na.rm = TRUE), "\n")
cat("FDR significant (q < 0.05):", sum(summary_table$sig_fdr, na.rm = TRUE), "\n")
cat("Results: results/primary/ | Plots: results/figures/\n")
cat("========================================\n")
```

---

## 9. STROBE-MR REPORTING CHECKLIST

Follow Skrivankova et al. 2021, BMJ (doi:10.1136/bmj.n2233). Every item must appear in manuscript methods or supplementary:

- [ ] Population: European ancestry — all datasets
- [ ] IV selection: p < 5x10^-8, r2 < 0.001, F-statistic > 10 per SNP
- [ ] Clumping: r2 < 0.001, 10,000 kb window
- [ ] Three MR assumptions stated and evaluated
- [ ] All 5 MR methods reported
- [ ] Cochran's Q heterogeneity reported
- [ ] MR-Egger intercept pleiotropy reported
- [ ] MR-PRESSO outlier test reported
- [ ] Bidirectional MR for all significant associations
- [ ] MVMR for outcomes with >=2 significant metals
- [ ] All GWAS IDs and sample sizes reported
- [ ] Multiple testing stated: Bonferroni (p < 0.05/36 = 0.00139) + BH-FDR
- [ ] Power limitation noted for CSF GWAS (n~3,146)
- [ ] R package versions reported (sessionInfo())

---

## 10. KEY LIMITATIONS (Address in Discussion)

| Limitation | How to Address |
|------------|----------------|
| CSF GWAS n~3,146 | Limits power; motivate replication with larger future cohorts |
| European ancestry only | Limits generalizability; state explicitly |
| Proxy GWAS for metals | Some metal GWAS use metabolomic proxies not direct metal measurement |
| Horizontal pleiotropy | Addressed by MR-PRESSO and MR-Egger intercept |
| Temporal ambiguity | MR estimates lifetime effect, not acute exposure |
| Potential sample overlap | Check for exposure-outcome GWAS overlap; flag if present |

---

## 11. TARGET JOURNALS (Ranked)

| Journal | IF | Rationale |
|---------|----|-----------|
| Molecular Neurodegeneration | ~15 | Top tier; use if copper/iron story is strong; connects to EB1A copper-AD hypothesis |
| Alzheimer's Research & Therapy | ~8 | Strong fit; publishes systematic MR + biomarker studies |
| Translational Psychiatry | ~6 | Published similar MR/AD biomarker designs |
| Neurobiology of Aging | ~4.5 | Good for mechanistic framing |
| Journal of Alzheimer's Disease | ~4 | Reliable fallback |

---

## 12. KEY REFERENCES

1. Liu X et al. (2019). Mineral Nutrition and the Risk of Chronic Diseases: A Mendelian Randomization Study. Nutrients, 11(2), 378.
2. Meng L et al. (2022). Are micronutrient levels and supplements causally associated with the risk of Alzheimer's disease? Food & Function, 13(12), 6665-6673.
3. Fu LL et al. (2025). The role of 25-OH vitamin D in Alzheimer's disease through MR and MRI. QJM, 118(1), 24.
4. Desai V & Kaler SG (2023). Metals in Alzheimer's Disease. Biomedicines, 11(4), 1161. [Key review explicitly calling for metals x CSF biomarker MR studies]
5. Deming Y et al. (2017). GWAS identifies four novel loci associated with Alzheimer's endophenotypes. Acta Neuropathologica, 133(5), 839-856. [Source: CSF GWAS]
6. Hibar DP et al. (2017). Novel genetic loci associated with hippocampal volume. Nature Communications, 8, 13624. [Source: Hippocampal volume GWAS]
7. Sinnott-Armstrong N et al. (2021). Genetics of 35 blood and urine biomarkers in the UK Biobank. Nature Genetics, 53, 185-194. [Source: Copper, zinc, selenium GWAS]
8. Skrivankova VW et al. (2021). STROBE-MR. BMJ, 375, n2233. [Reporting standard]
9. Larsson SC et al. (2025). Human Trace Elements, Gut Microbiota, and Alzheimer's Disease: Multistage MR Analysis. Food Science & Nutrition. [Still uses clinical AD as outcome — confirms our gap]

---

## 13. BEHAVIOR RULES

1. Never fabricate rsIDs, betas, p-values, or GWAS statistics. If OpenGWAS returns an error, report it clearly and stop.
2. Always show intermediate outputs — SNP counts, F-stats, Q-stats, harmonization logs — visible at every stage.
3. Stop at failed checkpoints. Ask Ashwin before proceeding past any flagged stage.
4. Report all 5 MR methods always. IVW is never reported alone.
5. Never suppress warnings to make results look cleaner.
6. Apply BH-FDR across all 36 IVW p-values — mandatory.
7. Separate statistical and biological significance explicitly for every result.
8. Novelty check before writing. Confirm with Ashwin before claiming any finding is novel.
9. Manuscript sentences must trace to validated entries in logs/session_log.md.
10. Verify all GWAS IDs at the start of every analysis session before touching code.

---

## 14. SESSION START CHECKLIST

Confirm with Ashwin at the start of every session:
- [ ] Which stage are we running today?
- [ ] Which pairs are in scope?
- [ ] Have all GWAS IDs been verified (especially Sinnott-Armstrong metals and cortical thickness)?
- [ ] Has the pilot run been completed (iron → CSF Abeta42, end to end) before full batch run?
- [ ] Are there newer CSF GWAS to check (Schott 2021, Palmqvist 2021)?
- [ ] Continuing from prior session or fresh start?

Do not begin analysis until all confirmed.

---

## 15. MANUSCRIPT MODE

When Ashwin switches to writing:
- Ground every sentence in a result from logs/session_log.md
- Use MR-appropriate hedged language: "MR evidence supports a potential causal effect of X on Y" — never "X causes Y"
- Flag any causal claim without HIGH confidence tier support
- Target: Molecular Neurodegeneration — concise, mechanistic, citation-dense
- Required sections: Abstract / Introduction / Methods / Results / Discussion
- Methods must include: instrument selection, harmonization, all 5 MR methods, all sensitivity analyses, FDR correction method, R package versions

---

*Last updated: March 27, 2026. GWAS IDs must be re-verified at analysis time. Update this file as the study evolves.*
