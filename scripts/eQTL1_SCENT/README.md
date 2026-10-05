# eQTL1: SCENT-style gene bootstrap

The original script `run_eqtl1_scent.R` reports the pooled OR only.  The new
script `run_eqtl1_bootstrap_both.R` reports both the pooled OR and the
SCENT-style mean of per-gene enrichments using the same 1,000 bootstrap draws.

For each dataset, broad cell type, and prediction method:

1. Use the method's own `significant_pairs.tsv.gz` calls and retain positive
   enhancer-gene effects.
2. Define causal eQTLs as matched OneK1K SuSiE variant-gene pairs with
   PIP >= 0.20 and distance >1 kb from the gene TSS.
3. For each gene, count causal and 1000G background SNPs inside versus outside
   the selected linked peaks, within the common 500 kb candidate window.
4. Compute a pooled one-sided Fisher test and Haldane-Anscombe odds ratio.
5. Sample eligible genes with replacement 1,000 times and recompute both
   statistics.  This is what “bootstrap different genes” means: genes, not
   cells or individual links, are the independent resampling units.

For the SCENT-style statistic, each eligible gene first gets the causal-variant
enrichment used in the SCENT paper:

```text
(causal variants in linked peaks / background variants in linked peaks)
---------------------------------------------------------------------
(all causal variants in cis / all background variants in cis)
```

The reported point estimate is the arithmetic mean across genes.  Genes with
zero common background variants inside the linked annotation have undefined
enrichment and are omitted from that mean, following the paper's formula;
they are not given a pseudocount-based enrichment.

Run all configured datasets with:

```bash
Rscript scripts/eQTL1_SCENT/run_eqtl1_scent.R
```

Run both summaries with exactly 1,000 bootstrap replicates:

```bash
Rscript scripts/eQTL1_SCENT/run_eqtl1_bootstrap_both.R
```

The output is written to `results/validation/eQTL/SCENT_both/`.  Each method
has `gene_bootstrap_pooled_or_and_scent_mean.tsv.gz`, containing one row per
bootstrap replicate, and `summary.tsv`, containing point estimates and 95% CIs.

To run one dataset/cell type, append the config path, dataset, and cell type,
for example `config/eqtl_validation.yml Dogma_CTRL Tnk`.
