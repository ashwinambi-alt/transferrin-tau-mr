# Transferrin and Tau Biomarkers: Multivariable Mendelian Randomization

**Paper:** Iron transport capacity causally reduces circulating and CSF tau biomarkers independent of serum iron: a multivariable Mendelian randomization study

**Author:** Ashwin Ambi
**Affiliation:** Enzo Life Sciences
**Preprint:** [TO BE ADDED]
**Journal:** Submitted to Journal of Neurology, Neurosurgery & Psychiatry

---

## Overview

Two-sample Mendelian randomization study testing whether genetically predicted transferrin protein concentration causally affects tau biomarkers in blood and cerebrospinal fluid (CSF), independent of serum iron. Primary finding: MVMR beta = -0.044, SE = 0.010, p = 2.29e-05.

---

## Repository Contents

| File | Description |
|------|-------------|
| `MASTER_PIPELINE.R` | Full TwoSampleMR analysis pipeline covering 36 primary analyses |
| `REPRODUCIBILITY_TEST.R` | 15-run stability validation confirming numerical reproducibility |
| `CLAUDE.md` | Pipeline documentation (v3.0) |
| `validation_audit_primary.csv` | Fresh API replication of primary results (March 2026) |
| `validation_audit_CSF.csv` | CSF subgroup replication results |
| `TableS1_SNP_instruments.csv` | Genetic instruments table (8 SNPs) |
| `TableS2_complete_gradient.csv` | Complete gradient across all outcome pairs tested |

---

## How to Run

1. Open `MASTER_PIPELINE.R` in RStudio
2. Install required packages if not already installed:
```r
install.packages(c("TwoSampleMR", "MendelianRandomization", "dplyr", "ggplot2", "MRPRESSO"))
```

3. Source the pipeline — all GWAS data are accessed live via the IEU OpenGWAS API. No local data files are required.
4. Run `REPRODUCIBILITY_TEST.R` to verify numerical stability across 15 seeds.

**R version:** 4.5.3 (2026-03-11 ucrt)
**TwoSampleMR version:** 0.7.1
**MendelianRandomization version:** 0.10.0

---

## Data Sources

All GWAS summary statistics accessed from:
- IEU OpenGWAS: https://opengwas.io/
- GWAS Catalog: https://www.ebi.ac.uk/gwas/

No individual-level data were used. Key accessions:
- Transferrin exposure: ieu-a-1052 (Benyamin 2014, n=23,986)
- Primary outcome: ebi-a-GCST90095138 (Sarnowski 2022, n=14,721)
- CSF outcomes: GCST90129600, GCST90134631, GCST90134632 (Jansen 2022)

---

## Citation

> Ambi A. Iron transport capacity causally reduces circulating and CSF tau biomarkers independent of serum iron: a multivariable Mendelian randomization study. [Journal], [Year]. [DOI TO BE ADDED]

---

## License

MIT License — code is freely available for reuse with attribution.
