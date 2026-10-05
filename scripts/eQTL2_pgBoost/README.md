# eQTL2: pgBoost/eQTLEnrichment-style ranking

This analysis evaluates continuous ranking separately from method-specific
significance calibration.

1. OneK1K variant-gene pairs with SuSiE PIP >= 0.50 are intersected with the
   common candidate peak universe.  A candidate SNP-peak-gene unit is positive
   only when the fine-mapped eGene equals the predicted target gene.
2. Each method's positive, signed score is ranked.  The primary EpiLink and GLM
   score is `beta`; scMultiMap uses `cor`; Signac uses `score`.  A parallel
   standardized-evidence analysis uses `zval`, `test_stat`, or `zscore`.
3. At fixed top-rank fractions, the pipeline reports recall, precision, and
   enrichment over the assessable base rate.  It also summarizes mean
   enrichment across a prespecified recall grid.
4. The method's native positive significant set is evaluated separately.  This
   is the result that can differ when EpiLink changes FDR calibration while the
   underlying fitted coefficients are unchanged.

EpiLink and its paired Poisson/NB GLM consume identical `all_pairs` estimates
in this project.  Therefore their beta-ranked curves are expected to overlap
exactly.  Their native-significance points can differ because EpiLink uses the
empirical-null FDR rule and the GLM baseline uses conventional BH selection.

Run all configured datasets with:

```bash
Rscript scripts/eQTL2_pgBoost/run_eqtl2_ranked.R
```
