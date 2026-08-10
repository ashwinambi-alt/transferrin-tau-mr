# Transferrin and Tau Biomarkers: Multivariable Mendelian Randomization

**Paper:** Genetically predicted transferrin level is associated with lower circulating total tau independent of serum iron: a Mendelian randomization study

**Author:** Ashwin Ambi
**Affiliation:** Enzo Life Sciences, Inc., Farmingdale, NY, USA
**Journal:** Journal of Alzheimer's Disease (major revision, 2026)
**Preprint:** [TO BE ADDED]

---

## Overview

Two-sample Mendelian randomization study testing whether genetically predicted serum transferrin is associated with tau biomarkers in blood and cerebrospinal fluid (CSF), independent of serum iron. The findings are framed as **associational / hypothesis-generating**, not confirmatory causal claims.

**Primary finding:** genetically predicted transferrin is inversely associated with circulating total-tau (IVW β = −0.039, SE = 0.013, p = 2.88×10⁻³, 8 SNPs); multivariable MR conditioning on serum iron does not attenuate the estimate (β = −0.044, SE = 0.009, **p = 1.4×10⁻³**, t-based; serum iron null, p = 0.80). The signal is reproduced by a cis-restricted TF-locus analysis (β = −0.060, p = 3.99×10⁻⁴) and supported by colocalization (conditional PP.H4/(PP.H3+PP.H4) = 96.3%). The univariable result does not survive Bonferroni correction across the exploratory outcome panel.

---

## Repository Contents

### Original submission
| File | Description |
|------|-------------|
| `MASTER_PIPELINE.R` | Full TwoSampleMR analysis pipeline (36 primary analyses) |
| `REPRODUCIBILITY_TEST.R` | 15-run stability validation |
| `validation_audit_primary.csv` / `validation_audit_CSF.csv` | Fresh API replication (March 2026) |
| `TableS1_SNP_instruments.csv` | Genetic instruments (8 SNPs) |
| `TableS2_complete_gradient.csv` | Complete gradient across all outcome pairs |

### Major revision (`revision_R1/`)
Reviewer-response analyses added for the JAD revision: cis-restricted TF-locus MR, Bayesian colocalization, leave-locus-out (random-effects IVW), per-instrument GWAS Catalog pleiotropy audit, full MVMR reporting (Sanderson–Windmeijer conditional F, Q_A, MVMR-Egger, 3-/4-way models), and transferrin-specific positive/negative controls. See `revision_R1/REVISION_NOTES.md` for the file index and the corrected MVMR p-value note (t-based 1.4×10⁻³ supersedes the earlier z-approximation 2.29×10⁻⁵).

---

## How to Run

1. Open `MASTER_PIPELINE.R` in RStudio; install packages:
```r
install.packages(c("TwoSampleMR", "MendelianRandomization", "MVMR", "coloc", "dplyr", "ggplot2", "MRPRESSO"))
```
2. Most GWAS data are accessed live via the IEU OpenGWAS API. Revision scripts take the project directory as their first argument.

**R version:** 4.5.3 / 4.6.1 · **TwoSampleMR:** 0.7.1 · **MendelianRandomization:** 0.10.0 · **MVMR** (t-based inference) · **coloc** (colocalization)

---

## Data Sources

All GWAS summary statistics from IEU OpenGWAS (https://opengwas.io/) and the GWAS Catalog (https://www.ebi.ac.uk/gwas/). No individual-level data. Key accessions:
- Transferrin exposure: ieu-a-1052 (Benyamin 2014, n = 23,986)
- Primary outcome: ebi-a-GCST90095138 (Sarnowski 2022, n = 14,721)
- CSF outcomes: Jansen 2022 (GCST90129600 etc.)

---

## Citation

> Ambi A. Genetically predicted transferrin level is associated with lower circulating total tau independent of serum iron: a Mendelian randomization study. Journal of Alzheimer's Disease (under revision), 2026. [DOI TO BE ADDED]

## License

MIT License — code is freely available for reuse with attribution.
