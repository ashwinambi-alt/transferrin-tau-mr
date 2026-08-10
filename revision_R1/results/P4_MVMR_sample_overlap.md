# Part 4 — MVMR sample overlap: methodological response + Bell add-on plan

Approach (b): address the exposure–exposure sample-overlap concern methodologically now,
with the non-overlapping Bell 2021 iron MVMR to be added as a confirmatory sensitivity once
the deCODE file is obtained (see end).

## The overlap, stated explicitly

| Pair | Source | Overlap fraction | Consequence |
|---|---|---|---|
| Transferrin ↔ serum iron (the two MVMR **exposures**) | both Benyamin 2014, n=23,986 | **~100% (same cohort)** | affects SNP-exposure covariance / SEs, **not** point-estimate consistency (see below) |
| Exposure(s) ↔ total-tau **outcome** | Benyamin 2014 vs Sarnowski 2022 | **partial** — the manuscript Limitations discloses a shared Rotterdam Study sub-cohort | small bias toward the observational association; already disclosed |

**Estimated exposure–exposure overlap fraction: ~1.0** (transferrin and serum iron are measured
in the identical Benyamin 2014 individuals).

## Why the existing 2-way MVMR is still interpretable

For MVMR-IVW, complete overlap **between the exposures** does not bias the causal point
estimates provided the **outcome** GWAS is independent of the exposures (Sanderson et al. 2019;
Burgess et al. 2016). Overlap between exposures enters through the SNP-exposure covariance and
therefore affects the conditional-F and standard-error calculation, not the consistency of the
IVW estimator. The partial exposure–outcome overlap (Rotterdam) is the residual concern and is
disclosed as a limitation; it biases toward the observational (non-null) association rather than
creating a spurious one.

## gencov (SNP-exposure covariance) — what was done and the caveat

- The existing 2-way MVMR (Analysis E, `MVMR::ivw_mvmr` / `strength_mvmr`) used **gencov = 0**
  (the package default when a covariance matrix is not supplied).
- **Justification / caveat:** with the two exposures from the *same* Benyamin individuals, the
  true SNP-exposure estimation errors are correlated, so gencov = 0 is an approximation. It does
  not bias the point estimate; it can mildly mis-estimate the conditional-F and SEs. We report
  this transparently. Estimating a non-zero gencov properly requires bivariate LD-score
  regression intercepts or individual-level data, which are not available for these summary sets.
- **The clean resolution is the Bell add-on:** with serum iron from Bell 2021 (deCODE/Iceland +
  Denmark + UK) and transferrin from Benyamin, the two exposures are **independent samples**, so
  **gencov = 0 is exactly correct**, not an approximation.

## Existing 2-way MVMR numbers (Benyamin transferrin + Benyamin iron → tau)

Instrument set: union of transferrin- and serum-iron-associated SNPs, **10 SNPs** after allele
alignment. From `E_MVMR_extended.csv` / `E_MVMR_extended_summary.txt`:

- Transferrin: β = −0.0437, SE = 0.0092, **p = 1.4×10⁻³** (t-based, df = 10−2 = 8), conditional F = **278.2**
- Serum iron:  β = −0.0057, SE = 0.0218, **p = 0.801** (t-based, df = 8), conditional F = **93.9**
- MVMR heterogeneity Q_A = 6.28, **p = 0.508**
- SNP-effect correlation between exposures **r = −0.611**
- Both conditional F ≫ 10 → the null serum-iron coefficient is **not** an artefact of weak
  conditional instruments.

## Bell 2021 confirmatory MVMR — to run when the file is available

- **Serum iron source:** Bell 2021, GWAS Catalog **GCST011367, n = 163,511** (deCODE:
  `https://download.decode.is/form/2021/ironserum.7z`). *NB: the spec's "n = 246,139" is Bell's
  ferritin; true Bell serum iron is n = 163,511.*
- **Design:** re-run the 2-way MVMR with transferrin (Benyamin ieu-a-1052) + serum iron (Bell
  GCST011367) → total-tau (Sarnowski). Report: instrument set/count after alignment,
  Sanderson–Windmeijer conditional F for both, transferrin β/SE/p (t-based, state df), serum-iron
  β/SE/p, Q_A and p, SNP-effect correlation r. gencov = 0 is exactly justified (independent
  exposure samples); exposure–exposure overlap fraction = 0.
- **Expectation / interpretation:** if transferrin's estimate holds with an independent iron GWAS,
  the finding is robust to the exposure-overlap concern. If it changes, report plainly.
- Script stub to add: `R2_PART4_bell_mvmr.R` (reads local Bell iron sumstats; not yet run —
  awaiting the deCODE download).

## Bottom line for the response letter
The 2-way MVMR conclusion — transferrin's association with tau is not attenuated by conditioning
on serum iron — rests on strong conditional instruments and is not biased by the exposure–exposure
overlap under the standard MVMR result. We disclose the gencov = 0 approximation and the partial
exposure–outcome (Rotterdam) overlap, and we add a fully non-overlapping Bell 2021 iron MVMR as a
confirmatory sensitivity.
