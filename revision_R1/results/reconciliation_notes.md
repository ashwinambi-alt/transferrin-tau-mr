# Part 1 — reconciliation notes

Generated 2026-08-07. Source of truth for the primary panel: results/primary/harmonised_transferrin_tau_2026-03-28.csv (exposure ieu-a-1052 = Benyamin 2014 transferrin n=23,986; outcome ebi-a-GCST90095138 = Sarnowski 2022 total-tau n=14,721).

## 1.1 Cis-TF instrument count — RESOLVED: two different analyses, zero shared SNPs
- Fig S9 (Analysis C, script R1_C_cis_TF_MR.R): DE NOVO extraction of transferrin GWAS SNPs at p<5e-8 within +/-500 kb of the TF gene body (chr3:133,464,073-133,494,388), clumped r2<0.001 in EUR. -> 2 SNPs {rs3811658, rs17376530}, IVW beta=-0.060, SE=0.017, p=3.99e-4. Identical at +/-1 Mb.
- Fig S6 (Analysis B, script R1_B_leave_locus_out.R) row 'Cis-TF only (3 SNPs)': NOT a de novo extraction. It is the subset of the ORIGINAL 8-SNP primary panel hand-labelled TF/TFRC-region = {rs8177240, rs9990333, rs7646473}, IVW beta=-0.044, p=0.006.
- They share NO SNPs and are built differently. rs8177240 does not survive C's clumping because rs3811658 (849 bp away, high LD) tags the same signal and is retained instead.
- ACTION: report Analysis C as THE cis-only MR. Relabel the Fig S6 row 'TF/TFRC-region subset of primary panel (3 SNPs)', present as sensitivity, not a second cis definition.

## 1.2 Excl. HLA + FADS (6 SNPs) — RESOLVED: fixed- vs random-effects IVW
- Point estimate identical: beta=-0.0398 across 6 SNPs {rs1495741, rs1800562, rs744653, rs7646473, rs8177240, rs9990333}.
- Fixed-effect IVW: SE=0.0133, p=0.0029 (this is what Fig S6 / base-R Analysis B reports ~0.003).
- Multiplicative random-effect IVW: SE=0.0145, p=0.0060 (this matches manuscript/Fig S1 ~0.006).
- Cause: residual heterogeneity (Q=5.87, Q_p=0.319, I2=15%) inflates the RE SE. Not a harmonization difference and not a mislabeled row.
- ACTION: standardize on ONE convention. The manuscript primary (all-8) uses TwoSampleMR default IVW = multiplicative random-effects; for all-8 there is no heterogeneity so FE=RE and both give 2.88e-3. For internal consistency, report the leave-locus-out table with RANDOM-EFFECTS IVW (the conservative, manuscript-consistent choice) -> excl_HLA_FADS p=0.006. Regenerate Fig S6 accordingly.

## 1.3 Serum-iron MVMR p — RESOLVED: use 0.801 (t-based), retire 0.817 (z-based)
- 0.801 = MVMR package, t-distribution inference, df = n_snps - n_exposures = 10 - 2 = 8. SAME basis as the transferrin row (p=1.4e-3). Source: E_MVMR_extended.csv (serum_iron pval=0.8006).
- 0.817 = MendelianRandomization::mr_mvivw z-based approximation (old Fig 3). 0.80 = 0.801 rounded to 2 dp.
- ACTION: report 0.801 (or 0.80) everywhere incl. Fig 3 panels A & B. It is not acceptable to show transferrin t-based (1.4e-3) alongside serum iron z-based (0.817) in the same figure.

## 1.4 rs8177240 Wald SE — RESOLVED: 0.017
- From harmonized data: Wald beta = beta.outcome/beta.exposure = -0.0213/0.4234 = -0.0503; SE = se.outcome/|beta.exposure| = 0.0070/0.4234 = 0.01653 -> 0.017 (3 dp). p=0.0023.
- Fig S9 shows SE=0.016 (hardcoded 0.0165 that floors to 0.016). Manuscript's 0.017 is correct. ACTION: fix Fig S9 label to 0.017.

## 1.5 Colocalization posteriors — RESOLVED: sum to 1.000
- PP.H0=0.000  PP.H1=0.536  PP.H2=0.000  PP.H3=0.017  PP.H4=0.447  (sum=1.000)
- Fig S5 currently omits H0 and H2. ACTION: show all five to 3 dp; note PP.H4=0.447 is below the conventional 0.8 threshold and the conditional H4/(H3+H4)=96.3% is reported because PP.H1=0.536 reflects limited tau-GWAS power, not evidence against colocalization.

## Cross-cutting flag for Part 2
- Individual cis Wald ratios: rs3811658 beta=-0.0574 p=0.0014 (F=1269); rs17376530 beta=-0.0856 p=0.1039 (F=130).
- The cis IVW is ~90% weighted on rs3811658 (which tags rs8177240); rs17376530 is directionally concordant but individually non-significant (p~0.10). LD r2 (Part 2, API) needed to quantify independence. Report plainly.
