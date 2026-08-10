# Table S7 — transferrin-specific controls (Part 3)

Generated 2026-08-08. Exposure = 8-SNP transferrin instrument set (Benyamin ieu-a-1052), same pipeline as primary. TwoSampleMR IVW + sensitivity.
TIBC positive control unavailable in OpenGWAS (Bell TIBC = deCODE GCST011368 only); near-collinear with transferrin, so omitted with note.

## Results (IVW headline per outcome)
- **Transferrin saturation** (positive control): IVW b=-0.4810, p=0.000518 [8 SNPs]. Expect NEGATIVE (higher transferrin lowers saturation at fixed iron). Overlap: SHARES Benyamin 2014 cohort with the transferrin exposure.
- **Height** (negative control): IVW b=-0.0196, p=0.305 [8 SNPs]. Expect null. Overlap: independent (GIANT).
- **Educational attainment** (negative control): IVW b=-0.0041, p=0.222 [7 SNPs]. Expect null (no plausible iron pathway). Overlap: independent (SSGAC).
- **Serum iron** (DAG / mediator check): IVW b=-0.0915, p=0.524 [8 SNPs]. Expect tests whether iron is downstream of transferrin (mediator). Overlap: SHARES Benyamin 2014 cohort with the transferrin exposure.

## Interpretation
- Positive control (transferrin saturation): a strong NEGATIVE estimate validates the transferrin instruments' biology (more transferrin -> more unsaturated binding sites -> lower saturation), i.e. the NTBI-buffering premise. NOTE: TSAT shares the Benyamin cohort, so this is a within-sample biological-consistency check, not out-of-sample.
- Negative controls (height, educational attainment): null estimates confirm the transferrin instruments are not pleiotropically associated with unrelated traits.
- DAG / mediator check (transferrin -> serum iron): a null or weak estimate supports treating serum iron as a PARALLEL exposure (conditioned on in MVMR), not a mediator on the transferrin->tau path; a strong estimate would raise a mediator concern. Steiger direction reported. Report the result honestly either way.
