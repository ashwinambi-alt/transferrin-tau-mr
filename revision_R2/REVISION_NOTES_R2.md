# Revision round 2 (ALZ-26-0616.R2) — analyses added at the second review

The second reviewer report raised four points. Three required new computation;
the fourth was a wording change in the Discussion. The manuscript was accepted
on 23 September 2026 after these were addressed.

Nothing in the primary analysis changed. The instrument sets, the univariable
IVW estimate, the MVMR point estimates and Figures 1–3 are all identical to the
first revision. What changed is how two supporting analyses are characterised.

---

## 1. Colocalization: unconditional posteriors and prior sensitivity

The reviewer objected that reporting the *conditional* PP.H4 of 96.3% overstated
the result, since the unconditional posterior for a shared causal variant (0.45)
is lower than that for a transferrin-only signal (0.54). He was right.

**Unconditional posteriors at the default priors** (p1 = p2 = 1×10⁻⁴, p12 = 1×10⁻⁵):

| | posterior |
|---|---|
| PP.H0 no signal | 0.000 |
| PP.H1 transferrin only | **0.536** |
| PP.H2 tau only | 0.000 |
| PP.H3 distinct causal variants | 0.017 |
| PP.H4 shared causal variant | **0.447** |

Conditional PP.H4/(PP.H3+PP.H4) = 0.963. The data are therefore at least as
compatible with a transferrin-only regional signal as with a shared one, and the
paper no longer claims otherwise.

**Prior sensitivity** (`coloc_prior_grid.csv`, `coloc_sensitivity_curve.csv`,
`results/FigS11` in the supplement). Interpolated in log10(p12) from the
100-point curve:

- PP.H4 overtakes PP.H1 at **p12 ≈ 1.2005×10⁻⁵**, about **20% above** the default
- PP.H4 crosses 0.5 at **p12 ≈ 1.2388×10⁻⁵** (23.9%) — a *different* point

> Both figures are interpolated, not read off the grid. Taking the first grid
> point at which a condition already holds overstates the threshold: on the
> coarse six-point grid the crossing appears to be at 5×10⁻⁵, and on the
> 100-point curve at 1.29×10⁻⁵. Neither is the crossing.

**Regional detail** (`coloc_regional_detail.csv`, `coloc_window_counts.csv`).
Window chr3:132,964,073–133,994,388 (GRCh37), ±500 kb of the *TF* gene body.
817 transferrin and 3151 total-tau variants before harmonisation; 696 and 3148
after QC (MAF > 0.005, finite SE); **696 shared and analysed**.

The regional lead variant in the tau GWAS is **rs3811658** (p = 1.44×10⁻³), not
rs8177240, which ranks third of 696 (p = 2.34×10⁻³). They lie 849 bp apart with
r² = 0.974, so no estimate changes, but the paper states it.

## 2. Cis instruments: per-variant Wald ratios and IVW weights

The cis-restricted analysis does not independently confirm the rs8177240 signal,
because one of its two instruments tags it.

| variant | Wald β (SE) | p | IVW weight | r² with rs8177240 |
|---|---|---|---|---|
| rs3811658 | −0.0574 (0.0180) | 0.0014 | **89.5%** | 0.974 |
| rs17376530 | −0.0856 (0.0526) | 0.104 | **10.5%** | 0.0013 |

Pooled cis IVW β = −0.0604, SE = 0.0171, p = 3.9892×10⁻⁴. rs17376530 — the only
genuinely independent TF-region variant — is directionally concordant but not
individually significant (95% CI −0.189 to 0.018).

The cis p-value is no longer described as "clearing Bonferroni". It falls below
the threshold, but the analysis was added post hoc at review and re-estimates the
same TF-locus signal.

> `TableS6_with_weights.csv` carries the weights. Recomputing the pooled estimate
> from the 3-dp values printed in Supplementary Table S6 gives p = 3.995×10⁻⁴;
> the harmonised data is authoritative.

## 3. MVMR between-exposure covariance

Both exposures come from Benyamin 2014, so the exposure-exposure sample overlap
is complete and the between-exposure covariance was set to zero. The reviewer
asked for that assumption to be estimated or tested.

Benyamin 2014 is a meta-analysis across eleven discovery cohorts, so no
sample-specific phenotypic correlation is available. A sensitivity analysis was
run instead across assumed phenotypic correlations ρ, converted to a covariance
matrix with `MVMR::phenocov_mvmr` (`mvmr_covariance_sensitivity.csv`):

| ρ | cond. F transferrin | cond. F serum iron | Q_A (7 df) | p |
|---|---|---|---|---|
| −0.8 | 1262.8 | 142.3 | 6.285 | 0.507 |
| −0.4 | 455.9 | 113.1 | 6.283 | 0.507 |
| **0 (assumed)** | **278.2** | **93.9** | **6.280** | **0.508** |
| +0.4 | 200.2 | 80.2 | 6.277 | 0.508 |

Conditional F never falls below 80 — eight times the conventional threshold of
10 — and Q_A stays null throughout.

The MVMR-IVW estimates are unchanged at every ρ, and this is algebraic rather
than empirical: `MVMR 0.4.6`'s `ivw_mvmr(r_input, gencov)` accepts the covariance
argument but never uses it. It warns when `gencov == 0` and then fits
`lm(betaYG ~ -1 + betaX1 + betaX2, weights = 1/seBetaYG^2)`. The covariance
enters only `strength_mvmr()` and `pleiotropy_mvmr()`.

Transferrin rises in iron deficiency, so the expected sign of ρ is negative — and
negative ρ *raises* the conditional F-statistics. Setting the covariance to zero
was the conservative choice.

---

## Files

```
revision_R2/
  results/
    coloc_prior_grid.csv              posteriors across the prior grid
    coloc_sensitivity_curve.csv       100-point p12 curve (source of the interpolated thresholds)
    coloc_regional_detail.csv         window, lead variants, LD, top SNP.PP.H4
    coloc_window_counts.csv           variants per trait before/after harmonisation
    TableS6_with_weights.csv          cis instruments with IVW weights
    TableS8_colocalization.csv        source for Supplementary Table S8
    mvmr_covariance_sensitivity.csv   conditional F and Q_A across rho
    R2_ANALYSES_console.log           full console output, including sessionInfo()
  scripts/
    R2_ANALYSES.R                     analyses 1-3; runs offline from revision_R1 outputs
    R2_A4_window_counts.R             the one step needing live OpenGWAS access
    figures/
      FIX_FIGS11_R2_style.R           Figure S11 (new)
      FIX_FIGS5_R2_label.R            Figure S5, relabelled for rs3811658
      FIX_FIGS9_R2_caption.R          Figure S9, caption corrected
```

`R2_ANALYSES.R` recomputes the colocalization posteriors from the per-SNP
approximate Bayes factors saved in `revision_R1/results/D_coloc_TF_snp_posteriors.csv`,
so it needs no network access; it reproduces the published posteriors to 5×10⁻¹⁴.
Scripts take `PROJECT_DIR` as the first argument.

**Versions:** R 4.6.1 · TwoSampleMR 0.7.9 · MendelianRandomization 0.10.0 ·
MVMR 0.4.6 · coloc 5.2.3 — unchanged from the first revision.
