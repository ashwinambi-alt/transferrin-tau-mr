# Revision R1 — reviewer-response analyses (Journal of Alzheimer's Disease)

Analyses added during the major revision, answering Reviewer 2's requests. All numbers below are canonical for the revised manuscript.

## Corrected MVMR p-value (transparency)
The original submission reported MVMR transferrin **p = 2.29×10⁻⁵**, from a z-based normal approximation (`MendelianRandomization::mr_mvivw`). The t-distribution inference appropriate to the residual degrees of freedom (`MVMR::ivw_mvmr`) gives an **identical point estimate**; the standard error changes only slightly between packages (0.010 vs 0.009), and the t-inference on the 8 residual degrees of freedom gives **p = 1.4×10⁻³**. The revised manuscript reports 1.4×10⁻³; the z-approximation value is superseded. The MVMR 95% confidence intervals use the same t-distribution (see "Final-round corrections" below).

## Key results
- Univariable IVW: β = −0.039, SE = 0.013, p = 2.88×10⁻³ (8 SNPs) — does not survive Bonferroni (0.00217).
- MVMR (2-way, adj. serum iron): transferrin β = −0.044, SE = 0.009, t-based 95% CI −0.065 to −0.022, p = 1.4×10⁻³; serum iron null, p = 0.80 (95% CI −0.056 to +0.045). Conditional F 278.2 / 93.9; Q_A p = 0.51; r = −0.611. Not presented as clearing Bonferroni: it re-estimates the same association conditional on serum iron.
- MVMR-Egger: intercept p = 0.24; transferrin slope β = −0.056, p = 1.65×10⁻⁴ (normal approximation).
- Cis-only TF (rs3811658 + rs17376530): β = −0.060, SE = 0.017, p = 3.99×10⁻⁴, clearing Bonferroni. LD: rs8177240↔rs3811658 r² = 0.974; rs17376530↔rs8177240 r² = 0.0013. The estimate rests principally on rs3811658, which tags the rs8177240 signal; rs17376530 is directionally concordant but not individually significant (p = 0.10). It is a transferrin-specific re-estimation, not an independent replication.
- rs8177240: 63.5% of IVW weight; Wald β = −0.050, SE = 0.017, p = 2.34×10⁻³.
- Colocalization: PP.H1 = 0.54, PP.H3 = 0.02, PP.H4 = 0.45; conditional PP.H4/(PP.H3+PP.H4) = 96.3%.
- Leave-locus-out (random-effects IVW): joint HLA+FADS+HFE+NAT2 exclusion β = −0.043, WM p = 0.005, RE-IVW p = 0.053; excl. rs8177240 alone β = −0.020, p = 0.358.
- Transferrin controls: TSAT (positive; shares Benyamin cohort) β = −0.481, p = 5.18×10⁻⁴; height p = 0.31, education p = 0.22 (null); serum iron IVW p = 0.52, WM +0.16 p < 10⁻⁴.
- APOE4 × transferrin interaction p = 0.095 (ns); CSF Aβ42 interaction p = 0.96.

## Scripts (`scripts/`)
- `R1_A_..H_*.R` — rs8177240 weight, leave-locus-out, cis-only TF MR, colocalization, extended/4-way MVMR, GWAS Catalog pleiotropy, MVMR-Egger/APOE4 polish.
- `R1_AUDIT.R` — full audit vs manuscript + source CSVs.
- `R2_PART1..4_*.R` — reconciliation of conflicts, per-SNP cis + LD, transferrin controls, MVMR sample-overlap.
- `R1_FIGURES.R` renders Fig S6. It also contains earlier versions of S5 and S7–S9. `FIX_FIG3_print.R` and `R2_PART3_FigS10.R` are earlier versions of Fig 3 and Fig S10, **superseded** by the scripts in `scripts/figures/`.

### Final figures (`scripts/figures/`) — the scripts that rendered the submitted figures
| Figure | Script | Inputs |
|---|---|---|
| Fig 1, Fig S3 | `FIX_FIG1_FIGS3_exact.R` | full-precision values written into the script, with their sources listed in its header |
| Fig 2, Fig S1, Fig S2 | `FIX_FIG2_FIGS1_FIGS2.R` | `data/FINAL_analysis_results.RData` (harmonised data; not versioned, `*.RData` is git-ignored) |
| Fig 3 | `FIX_FIG3_tCI.R` | values written into the script (from `results/E_MVMR_extended.csv`) |
| Fig S4 | `FIX_FIGS4_caption.R` | values written into the script (CSF subgroup results) |
| Fig S5, Fig S7 | `FIX_FIGS7_tCI.R` | `D_coloc_TF_snp_posteriors.csv`, `E_MVMR_extended.csv`, `F_MVMR_4way.csv` |
| Fig S6 | `../R1_FIGURES.R` | `B_leave_locus_out.csv` |
| Fig S8 | `FIX_FIGS8_keysize.R` | none (causal diagram) |
| Fig S9 | `FIX_FIGS9_notation.R` | `C_cis_TF_mr_500kb.csv` |
| Fig S10 | `FIX_FIGS10_notation.R` | `TableS7_transferrin_controls.csv` |

The figure scripts are committed exactly as run. They expect a project directory in which the files in `revision_R1/results/` sit under `results/reviewer_response/`, and they take that directory (or, for Fig 1/3/S4, an output directory) as their first argument. `FIX_FIGS10_notation.R` hard-codes the author's local project path and R library path in lines 12–15. Edit those lines before running it elsewhere.

Scripts take `PROJECT_DIR` as the first argument. Sample-overlap note: the between-exposure genetic covariance is set to zero (enters conditional F and Q_A, not the IVW point estimate); exposure–outcome overlap (Benyamin/Sarnowski share the Rotterdam Study) is disclosed in the manuscript Limitations.

## Final-round corrections (September 2026)
- **MVMR confidence intervals are t-based**, matching their t-based p-values: 2-way models use 8 residual df, and the 3-way model uses 28 df (31 SNPs − 3). This applies to Fig 3, Fig S7 and Table S2. Univariable IVW, weighted-median and Wald intervals stay at β ± 1.96·SE, because TwoSampleMR computes those p-values with the normal distribution.
- **Fig 1 is rendered from full-precision estimates.** The earlier script had back-calculated SEs from Table S2's rounded CIs. Its "Entorhinal thickness (bilateral)" row, a plain average of the left and right results, is replaced by separate L and R rows.
- **Table S2 CIs were recomputed** from full-precision β ± 1.96·SE (12 rows changed by 0.001–0.003).
- **`TableS2_complete_gradient.csv` (repository root) corrected:**
  - MVMR row: SE 0.009 and p 1.41e-03 (t-based), with serum iron p = 0.801
  - SEs for right entorhinal thickness (0.031 → 0.033), hippocampal volume (0.660 → 9.311), multiple sclerosis (0.028 → 0.0004) and depression (0.029 → 0.0009)

## Outputs (`results/`)
Editable table data (Tables S3–S7), reconciliation notes, LD matrix, colocalization posteriors, and per-analysis summaries.
