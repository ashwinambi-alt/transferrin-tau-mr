# Audit Report — reviewer-response analyses
Run: 2026-07-01 23:17:33.606244

Summary: 63 PASS, 0 FAIL, 0 WARN, 8 INFO

## Section A — Reproduce manuscript headline numbers

| Check | Expected | Observed | Status |
|---|---|---|---|
| Primary IVW β (manuscript = −0.039) | -0.039 | -0.03927 | PASS |
| Primary IVW SE (manuscript = 0.013) | 0.013 | 0.01318 | PASS |
| Primary IVW p (manuscript = 2.88e-3) | 0.00288 | 0.00288 | PASS |
| Primary IVW nsnp (manuscript = 8) | 8 | 8 | PASS |
| MVMR transferrin β (manuscript = −0.044) | -0.044 | -0.04365 | PASS |
| MVMR transferrin p (manuscript = 2.29e-5) | 2.29e-05 | 2.287e-05 | PASS |
| MVMR serum iron p (manuscript = 0.817) | 0.817 | 0.8166 | PASS |
| Egger intercept p (manuscript = 0.202) | 0.202 | 0.202 | PASS |
| Cochran Q p (manuscript = 0.544) | 0.544 | 0.544 | PASS |
| I² (manuscript = 0%) | 0 | 0 | PASS |
| MR-PRESSO global p (manuscript = 0.57) | 0.57 | 0.578 | PASS |
| Steiger correct SNPs (manuscript = 8/8) | 8/8 | 8/8 | PASS |
| LOO removals losing significance (manuscript = 1 (rs8177240)) | 1 | 1 | PASS |
| Weighted median p (manuscript = 0.001) | 0.00143 | 0.001534 | PASS |
| Positive control iron→Hb β (manuscript = +0.704) | 0.704 | 0.7044 | PASS |
| Negative control iron→height p (manuscript = 0.191) | 0.191 | 0.191 | PASS |

## Section B — GWAS ID staleness

| Check | Expected | Observed | Status |
|---|---|---|---|
| ieu-a-1049 trait | Benyamin 2014 serum iron | Iron | INFO |
| ieu-a-1049 sample size (manuscript = 23986) | 23986 | 23986 | PASS |
| ieu-a-1050 trait | Benyamin 2014 ferritin | Ferritin | INFO |
| ieu-a-1050 sample size (manuscript = 23986) | 23986 | 23986 | PASS |
| ieu-a-1051 trait | Benyamin 2014 TSAT | Transferrin Saturation | INFO |
| ieu-a-1051 sample size (manuscript = 23986) | 23986 | 23986 | PASS |
| ieu-a-1052 trait | Benyamin 2014 transferrin | Transferrin | INFO |
| ieu-a-1052 sample size (manuscript = 23986) | 23986 | 23986 | PASS |
| ieu-b-5115 trait | Bell 2021 ferritin (large) | ferritin | INFO |
| ieu-b-5115 sample size (manuscript = 246139) | 246139 | 246139 | PASS |
| ebi-a-GCST90095138 trait | Sarnowski 2022 total tau | Circulating levels of total-tau | INFO |
| ebi-a-GCST90095138 sample size (manuscript = 14721) | 14721 | 14721 | PASS |
| ebi-a-GCST90028995 trait | Vuckovic 2020 haemoglobin (pos ctrl) | Mean corpuscular hemoglobin | INFO |
| ebi-a-GCST90029008 trait | Vuckovic 2020 (neg ctrl) | Height | INFO |

## Section C — Cross-check summary values against source CSVs

| Check | Expected | Observed | Status |
|---|---|---|---|
| A_rs8177240_contribution.csv: rs8177240 weight_pct (report = 63.5) | 63.5 | 63.53 | PASS |
| A_rs8177240_contribution.csv: rs1800562 weight_pct (report = 27.15) | 27.15 | 27.15 | PASS |
| A_rs8177240_contribution.csv: sum of weight_pct == 100 | 100 | 100.0047 | PASS |
| B_leave_locus_out.csv: joint HLA+FADS+HFE+NAT2 ivw_beta (report −0.043) | -0.043 | -0.04278 | PASS |
| B_leave_locus_out.csv: joint ivw_pval (report 0.0073) | 0.0073 | 0.007271 | PASS |
| B_leave_locus_out.csv: cis_TF_3snps ivw_beta (report −0.044) | -0.044 | -0.04395 | PASS |
| B_leave_locus_out.csv: cis_TF_3snps wm_pval (report 0.0037) | 0.0037 | 0.003727 | PASS |
| C_cis_TF_mr_500kb.csv: cis-TF IVW β (report −0.060) | -0.06 | -0.06039 | PASS |
| C_cis_TF_mr_500kb.csv: cis-TF IVW p (report 3.99e-4) | 0.000399 | 0.0003989 | PASS |
| C_cis_TF_mr_500kb.csv: cis-TF nsnp (report 2) | 2 | 2 | PASS |
| D_coloc_TF_summary.csv: PP.H4 (report 0.447) | 0.447 | 0.447 | PASS |
| D_coloc_TF_summary.csv: PP.H3 (report 0.017) | 0.017 | 0.0171 | PASS |
| D_coloc_TF_summary.csv: PP.H1 (report 0.536) | 0.536 | 0.536 | PASS |
| D_coloc_TF_summary.csv: sum PP ≈ 1 | 1 | 1 | PASS |
| E_MVMR_extended.csv: MVMR transferrin β (report −0.044) | -0.044 | -0.04365 | PASS |
| E_MVMR_extended.csv: transferrin cond-F (report 278.2) | 278.2 | 278.2 | PASS |
| E_MVMR_extended.csv: serum iron cond-F (report 93.9) | 93.9 | 93.88 | PASS |
| E_MVMR_extended.csv: serum iron p (report 0.80) | 0.8 | 0.8006 | PASS |
| F_MVMR_4way.csv: 3-way transferrin β (report −0.046) | -0.046 | -0.04554 | PASS |
| F_MVMR_4way.csv: 3-way transferrin p (report 0.014) | 0.014 | 0.01371 | PASS |
| F_MVMR_4way.csv: 3-way transferrin cond-F (report 22.7) | 22.7 | 22.75 | PASS |
| F_MVMR_4way.csv: 3-way ferritin p (report ~0.82) | 0.824 | 0.824 | PASS |
| F_MVMR_4way.csv: 3-way TSAT p (report ~0.83) | 0.833 | 0.833 | PASS |
| G_gwas_catalog_summary.csv: rs8177240 n_iron (report 5) | 5 | 5 | PASS |
| G_gwas_catalog_summary.csv: rs1800562 n_hits (report 183) | 183 | 183 | PASS |
| G_gwas_catalog_summary.csv: rs8177240 concern (report LOW) | LOW | LOW — few non-iron hits | PASS |
| H_MVMR_Egger.csv: MVMR-Egger intercept p (report 0.244) | 0.244 | 0.244 | PASS |
| H_MVMR_Egger.csv: MVMR-Egger transferrin slope β (report −0.056) | -0.056 | -0.05624 | PASS |
| H_MVMR_Egger.csv: transferrin slope p (report 1.65e-4) | 0.000165 | 0.0001653 | PASS |
| H_interaction_tests.csv: APOE4 p-tau interaction p_Q (report 0.095) | 0.0952 | 0.0952 | PASS |
| H_interaction_tests.csv: Aβ42 APOE4 interaction p_Q (report 0.96) | 0.96 | 0.96 | PASS |

## Section D — Independent from-scratch recompute

| Check | Expected | Observed | Status |
|---|---|---|---|
| harmonised CSV row count | 8 | 8 | PASS |
| IVW β from independent recompute (should = −0.039) | -0.039 | -0.03927 | PASS |
| IVW p from independent recompute (should = 2.88e-3) | 0.00288 | 0.00288 | PASS |
| rs8177240 IVW weight % recompute (should = 63.5) | 63.5 | 63.53 | PASS |
| IVW β excl rs8177240 recompute (should = −0.020) | -0.02 | -0.02005 | PASS |
| IVW p excl rs8177240 recompute (should = 0.358) | 0.358 | 0.358 | PASS |
| harmonise action=2: effect_allele.exposure == effect_allele.outcome per SNP | TRUE | TRUE | PASS |
| beta.exposure sign distribution | 3 negative (rs1800562, rs7646473, rs9990333) | 3 negative (rs1800562, rs7646473, rs9990333) | PASS |
| cached excl-rs8177240 β matches manuscript −0.020 | -0.02 | -0.02005 | PASS |
| cached excl-rs8177240 p matches manuscript 0.358 | 0.358 | 0.358 | PASS |
