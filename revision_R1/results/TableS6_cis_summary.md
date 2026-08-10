# Table S6 — per-SNP cis instrument detail

Generated 2026-08-07. TF gene body GRCh37 chr3:133464073-133494388. LD from 1000G EUR (OpenGWAS).

## Individual Wald ratios
- rs8177240 (primary panel lead): beta=-0.0503, SE=0.0165, p=0.0023, F=1566, 0 bp from TF (inside gene body if 0).
- rs3811658 (cis instrument 1, tags rs8177240): beta=-0.0574, SE=0.0180, p=0.0014, F=1269, 0 bp from TF.
- rs17376530 (cis instrument 2, candidate independent): beta=-0.0856, SE=0.0526, p=0.1039, F=130, 60767 bp from TF.

## LD (1000G EUR)
- r2(rs8177240,rs3811658) = 0.9738
- r2(rs8177240,rs17376530) = 0.0013
- r2(rs3811658,rs17376530) = 0.0006

## Verdict on independence (bar: r2 < 0.01)
- rs3811658 vs rs8177240: r2 = 0.9738  -> NOT independent (expected high: same signal, 849 bp apart).
- rs17376530 vs rs8177240: r2 = 0.0013  -> INDEPENDENT of rs8177240.

## Plain-language conclusion
rs17376530 is INDEPENDENT of rs8177240 and supports the SAME (protective) direction (beta=-0.086), but is NOT individually significant (p=0.10). The cis-only IVW (beta=-0.060, p=3.99e-4) is ~90% weighted on rs3811658, which tags rs8177240. So the cis analysis is NOT simply rs8177240 relabelled (rs17376530 is a distinct, concordant variant), but neither is it two independently-significant instruments. Report this honestly: the cis result corroborates direction via an independent variant while its significance still leans on the rs8177240-tagging signal.

## Data-version note (rs8177240 / rs3811658 rsID labelling) — RESOLVED
- Harmonized/primary rs8177240 beta.exposure = 0.4234 (F=1566); the live ieu-a-1052 API returns 0.3883 for the rs8177240 label.
- RESOLVED against Ensembl GRCh37 (authoritative): rs8177240 = chr3:133,477,701, rs3811658 = chr3:133,476,852. The CACHED/manuscript data match Ensembl; the live OpenGWAS ieu-a-1052 API SWAPS these two adjacent rsID labels. Manuscript is correct; no estimate changes. See `P2_rsID_swap_investigation.md`.
- Reproducibility: map instruments by position, not live rsID, when using ieu-a-1052.
