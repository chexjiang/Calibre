# EpiLink eQTL validation

## Scope

The benchmark covers PBMC, BMMC, DOGMA-seq, arthritis tissue and SHARE-seq
populations listed in `config/eqtl_validation.yml`.  The current `pbmc` and
`bmmc` entries use the replacement QC'd inputs; the former `_new` suffix has
been retired.  Arthritis `T`, `B` and `Myeloid` use the OneK1K `Tnk`, `Bcell`
and `Myeloid` references, respectively.  SHARE-seq `LCL` uses the closest
available `Bcell` reference as an approximate match, since LCL is not a
primary OneK1K cell state.  Missing
method outputs are skipped and reported; they are never converted to negative
predictions.

## Shared reference preparation

1. Download eQTL Catalogue release-7 OneK1K SuSiE credible sets for B
   intermediate/memory/naive, CD14/CD16 monocytes, CD4/CD8 naive T cells, and
   NK cells.
2. Map Ensembl gene IDs to HGNC symbols with `EnsDb.Hsapiens.v86`.
3. Parse GRCh38 variant coordinates and collapse duplicate variant-gene pairs
   within each broad cell group by maximum PIP.
4. Use pgBoost's released 1000 Genomes MAC>5 SNP BED as the common-variant
   background.  Convert BED coordinates to one-based GRCh38 coordinates.
5. Exclude fine-mapped variants within 1 kb of the target TSS and use the same
   500 kb cis window as EpiLink candidate-pair generation.

Raw-file MD5 hashes and processed row counts are under
`reference/eQTL/shared/processed/`.

## eQTL1: SCENT-style enrichment

The prediction annotation for each method is its positive, native significant
peak-gene set.  For each gene, the pipeline counts:

|                    | In linked peaks | Outside linked peaks |
|--------------------|-----------------|----------------------|
| SuSiE PIP >= 0.20  | causal hit      | causal non-hit       |
| 1000G MAC>5 SNP    | common hit      | common non-hit       |

Counts are pooled for a one-sided Fisher test.  Genes with at least one causal
variant are sampled with replacement 1,000 times; the pooled odds ratio is
recalculated in each bootstrap.  Thus, “bootstrap genes” means resampling genes
as independent units, not resampling cells, peaks, or eQTL variants.  Empty
prediction sets have undefined OR and are reported as NA.

## eQTL2: pgBoost/eQTLEnrichment-style ranking

SuSiE PIP >= 0.50 variant-gene pairs are intersected with candidate peaks.  A
SNP-peak-gene unit is positive only if the predicted target equals the eGene.
At fixed top-rank fractions, the pipeline reports recall, precision, and
precision enrichment relative to the assessable base rate.

The primary continuous ranks are:

| Method | Primary rank | Evidence sensitivity rank |
|--------|--------------|---------------------------|
| EpiLink Poisson/NB | `abs(beta)` | `abs(zval)` |
| Poisson/NB GLM | `abs(beta)` | `abs(zval)` |
| scMultiMap | signed `cor` | signed `test_stat` |
| Signac | signed `score` | signed `zscore` |

EpiLink and its paired GLM use absolute beta/z-value ranks.  This matches the
two-sided EpiLink significance calibration and avoids discarding negative
links, since the OneK1K causal-link labels do not encode effect direction.
Ranks, rather than raw values, make score scales comparable across methods.

EpiLink and its paired GLM use the same fitted `all_pairs` coefficients in this
project.  Their beta- and z-ranked curves must therefore coincide.  The separate
native-significance benchmark can differ because EpiLink and the GLM select
different links through empirical-null versus conventional FDR calibration.

The two validation questions are intentionally separate.  eQTL1 evaluates the
method's native significant peak-gene set as a static enrichment test, whereas
eQTL2 evaluates continuous-score prioritization over recall.  Both use the same
OneK1K/SuSiE source, but eQTL1 uses PIP >= 0.20 plus a 1000G common-variant
background, while eQTL2 uses PIP >= 0.50 assessable SNP-peak-gene units and
their empirical base rate.  Signac corrected and uncorrected calls are compared
in the native-significance analyses: corrected uses `significant_pairs_corrected`
and uncorrected uses `significant_pairs_uncorrected`.  Their continuous `score`
ranking is identical, so only one Signac curve is shown in the eQTL2 ranking
plot.

## Commands

```bash
bash scripts/eQTL_shared/download_references.sh
Rscript scripts/eQTL_shared/prepare_reference.R
Rscript scripts/eQTL1_SCENT/run_eqtl1_scent.R
Rscript scripts/eQTL2_pgBoost/run_eqtl2_ranked.R
Rscript scripts/eQTL_shared/plot_eqtl_validation.R
Rscript scripts/eQTL1_SCENT/plot_eqtl1_bootstrap_both.R
Rscript scripts/eQTL_shared/plot_eqtl_dataset_level.R
```

The eQTL1 panel plots display Arthritis `B` and SHARE-seq `LCL` under the
common `Bcell` label.  OR is shown on a log10 axis; the SCENT-style mean uses a
pseudo-log axis so zero lower bootstrap limits remain visible.  The dataset-
level plots pool all configured cell types within each dataset: eQTL1 pools
gene-by-cell-type count units, while eQTL2 pools TP, selected units, and
positive units at each ranking fraction.
SHARE-seq/LCL is retained in eQTL1; its eQTL2 panel is omitted when no
assessable OneK1K PIP >= 0.50 units are available.

## Main outputs

The cleaned deliverable contains figures only:

- `results/validation/eQTL/OneK1K/eQTL1/`
- `results/validation/eQTL/OneK1K/eQTL2/`
- `results/validation/eQTL/OneK1K/figures/`

Interpret very high SCENT-style ORs with their bootstrap intervals and raw
counts: sparse fine-mapped hits can generate unstable ratios.  The ranking and
native-significance analyses answer different questions and should both be
shown rather than selecting whichever favors one method.
