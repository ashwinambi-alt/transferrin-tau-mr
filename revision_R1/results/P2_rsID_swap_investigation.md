# Forensic note: rs8177240 / rs3811658 rsID labelling at the TF locus

**Question raised (Part 2):** the live OpenGWAS API returned different exposure betas for
rs8177240 (0.388) and rs3811658 (0.423) than the cached harmonized files (rs8177240 = 0.423).
Which is correct?

**Investigation (2026-08-08).** Pulled full detail for both variants from OpenGWAS
(`ieu-a-1052` exposure and `ebi-a-GCST90095138` outcome) and adjudicated against the
authoritative dbSNP/Ensembl GRCh37 mapping (Ensembl REST, independent of OpenGWAS).

| rsID | Ensembl GRCh37 pos (authoritative) | Ensembl alleles | Cached pos / β_exp | Live OpenGWAS pos / β_exp |
|---|---|---|---|---|
| rs8177240 | **133,477,701** | T/C/G | 133,477,701 / 0.4234 ✓ | 133,476,852 / 0.3883 ✗ |
| rs3811658 | **133,476,852** | C/G/T | 133,476,852 / 0.3883 ✓ | 133,477,701 / 0.4234 ✗ |

**Conclusion — the cached/manuscript data is CORRECT; the live OpenGWAS `ieu-a-1052`
API swaps the rs8177240 and rs3811658 rsID labels between these two adjacent positions.**

Three independent confirmations that the cache is right:
1. Cached positions match Ensembl GRCh37 exactly for both variants.
2. Cached alleles are consistent with Ensembl (rs8177240 G/T ⊂ T/C/G; rs3811658 T/C ⊂ C/G/T).
3. The cached inter-variant distance is 133,477,701 − 133,476,852 = **849 bp**, matching the
   manuscript's statement that "rs3811658 sits 849 bp from rs8177240."

The two OpenGWAS datasets (exposure `ieu-a-1052`, outcome `ebi-a-GCST90095138`) agree with
each other but disagree with Ensembl — i.e. they share the same rsID-remap quirk at this
locus (the two variants are r²=0.97, 849 bp apart, and were evidently cross-assigned during
OpenGWAS's dbSNP harmonisation). Ensembl is authoritative and overrides.

## Impact
- **None on any estimate.** Every MR result in the manuscript uses the cached (correct)
  data. rs8177240 Wald: β = −0.0213 / 0.4234 = −0.050, SE = 0.007 / 0.4234 = 0.01653 → **0.017**
  (conflict 1.4 answer unchanged).
- rsIDs, positions, and the "849 bp" claim in the manuscript are all correct as written.
- The Part 2 LD result is unaffected: r² is looked up by rsID against 1000G and the
  rs8177240↔rs3811658 pair r²=0.974 / rs17376530↔rs8177240 r²=0.001 conclusions hold either way.

## Recommended action for the paper (reproducibility safeguard)
Add one sentence to the Methods/Data-availability or the code README:

> "Note: the OpenGWAS record `ieu-a-1052` (Benyamin 2014 transferrin) cross-labels the rsIDs
> rs8177240 (chr3:133,477,701) and rs3811658 (chr3:133,476,852) relative to dbSNP/Ensembl
> GRCh37. Our pipeline maps instruments by genomic position (Ensembl-authoritative); analysts
> reproducing this work via live rsID lookups on `ieu-a-1052` should map by position to avoid
> the swap."

This pre-empts a replicator or reviewer who re-runs against the live API and sees the labels
apparently disagree with the manuscript.
