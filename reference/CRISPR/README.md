# CRISPR-FlowFISH reference

These files are the validation references copied from
`incoming/validation/CRISPR`. The active FlowFISH implementation reads
`41586_2021_3446_MOESM7_ESM.txt` directly, retains its `CellType` column, and
maps the experiments to `Bcell`, `T`, `Myeloid`, and `K562` before lifting the
hg19 intervals to hg38 when needed. The pooled `crispr_hg38_positive.tsv` and
`crispr_hg38_negative.tsv` files are retained as legacy/archival files; they
are not required by `scripts/run_crispr_validation.R`.

Additional compatibility labels are `CD4T`/`Tnk` (Jurkat), `B` (B-cell
lines), and `LCL` (GM12878); `K562` remains a separate reference rather than
being merged into `Myeloid`.

Both EpiLink finite-cell correction settings are written separately under
`results/validation/CRISPR/<dataset>_enrichment_F.tsv` and
`results/validation/CRISPR/<dataset>_enrichment_T.tsv`. Each
correction directory contains one `enrichment_metrics.tsv`; its `celltype`
column identifies matched dataset cell types plus `All`, which is formed from
the union of the matched cell types only.
