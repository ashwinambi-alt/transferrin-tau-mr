# Revision R1 — reviewer-response analyses (Journal of Alzheimer's Disease)

Analyses added during the major revision, answering Reviewer 2's requests. All numbers below are canonical for the revised manuscript.

## Corrected MVMR p-value (transparency)
The original submission reported MVMR transferrin **p = 2.29×10⁻⁵**, from a z-based normal approximation (`MendelianRandomization::mr_mvivw`). The t-distribution inference appropriate to the residual degrees of freedom (`MVMR::ivw_mvmr`) gives the **same point estimate and SE** but **p = 1.4×10⁻³**. The revised manuscript reports 1.4×10⁻³; the z-approximation value is superseded.

## Key results
- Univariable IVW: β = −0.039, SE = 0.013, p = 2.88×10⁻³ (8 SNPs) — does not survive Bonferroni (0.00217).
- MVMR (2-way, adj. serum iron): transferrin β = −0.044, SE = 0.009, p = 1.4×10⁻³; serum iron null p = 0.80. Conditional F 278.2 / 93.9; Q_A p = 0.51; r = −0.611.
- MVMR-Egger: intercept p = 0.24; transferrin slope β = −0.056, p = 1.65×10⁻⁴.
- Cis-only TF (rs3811658 + rs17376530): β = −0.060, SE = 0.017, p = 3.99×10⁻⁴. LD: rs8177240↔rs3811658 r² = 0.974; rs17376530↔rs8177240 r² = 0.0013.
- rs8177240: 63.5% of IVW weight; Wald β = −0.050, SE = 0.017, p = 2.34×10⁻³.
- Colocalization: PP.H1 = 0.54, PP.H3 = 0.02, PP.H4 = 0.45; conditional PP.H4/(PP.H3+PP.H4) = 96.3%.
- Leave-locus-out (random-effects IVW): joint HLA+FADS+HFE+NAT2 exclusion β = −0.043, WM p = 0.005, RE-IVW p = 0.053; excl. rs8177240 alone β = −0.020, p = 0.358.
- Transferrin controls: TSAT (positive; shares Benyamin cohort) β = −0.481, p = 5.18×10⁻⁴; height p = 0.31, education p = 0.22 (null); serum iron IVW p = 0.52, WM +0.16 p < 10⁻⁴.
- APOE4 × transferrin interaction p = 0.095 (ns); CSF Aβ42 interaction p = 0.96.

## Scripts (`scripts/`)
- `R1_A_..H_*.R` — rs8177240 weight, leave-locus-out, cis-only TF MR, colocalization, extended/4-way MVMR, GWAS Catalog pleiotropy, MVMR-Egger/APOE4 polish.
- `R1_AUDIT.R` — full audit vs manuscript + source CSVs.
- `R1_FIGURES.R` — supplementary figures S5–S9. `FIX_FIG3_print.R` — main Fig 3. `R2_PART3_FigS10.R` — Fig S10 (controls).
- `R2_PART1..4_*.R` — reconciliation of conflicts, per-SNP cis + LD, transferrin controls, MVMR sample-overlap.

Scripts take `PROJECT_DIR` as the first argument. Sample-overlap note: the between-exposure genetic covariance is set to zero (enters conditional F and Q_A, not the IVW point estimate); exposure–outcome overlap (Benyamin/Sarnowski share the Rotterdam Study) is disclosed in the manuscript Limitations.

## Outputs (`results/`)
Editable table data (Tables S3–S7), reconciliation notes, LD matrix, colocalization posteriors, and per-analysis summaries.
