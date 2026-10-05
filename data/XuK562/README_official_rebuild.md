# Xu K562 official reconstruction

This directory is rebuilt from the official ISSAAC-seq/10x K562 runs rather than the previous partial h5mu object.

- RNA: ENA `ERR9847050`
- ATAC: ENA `ERR9847049`
- Reference workflow: `https://github.com/dbrg77/ISSAAC-seq`, `other_methods_preprocessing/10x_multiome/K562`
- QC thresholds (Dorans et al. Supplementary Table 14): RNA counts > 1000, ATAC counts > 200, nucleosome signal < 2, TSS enrichment > 2.
- 10x ARC uses different GEX and ATAC barcode sequences. `run_official_pipeline.sh` translates ATAC barcodes to GEX barcodes by matching the corresponding line in `atac_737K-arc-v1_rc.txt` and `gex_737K-arc-v1.txt` before pairing cells.

Run `run_official_pipeline.sh` after the SRA files are available. It runs STARsolo for RNA, Chromap for ATAC, barcode translation, and genome-wide MACS2 peak calling from an evenly spaced fragment subsample (MACS2's paired-end parser requires a seekable input). Then run `assemble_official_xuk562.R` to pair RNA/ATAC cells, calculate nucleosome/TSS metrics, apply the thresholds, and write the EpiLink RDS files. Finally run `data_preprocessing.R` to create `XuK562.gene_peak_df.rds`.

The current rebuild produced 11,861 RNA/ATAC barcode pairs before QC and 6,949 cells after QC. The files are `XuK562.rna.rds`, `XuK562.atac.rds`, `XuK562.meta.rds`, `XuK562.qc_summary.tsv`, and `XuK562.gene_peak_df.rds`.
