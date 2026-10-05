# Satpathy K562 multiome benchmark data

This directory contains the Satpathy K562 subset of the public IGVF multiplexed
multiome analysis set `IGVFDS2018URDY`. The released analysis set contains a
MULTI-seq mixture of Jurkat, K562, THP-1 and TeloHAEC cells, with an unfiltered
RNA HDF5 matrix, cell-hashing counts, GRCh38 ATAC fragments and a GRCh38 peak
set. The official files and metadata are available from the
[IGVF data portal](https://data.igvf.org/) (accession
`IGVFDS2018URDY`).

## Processing

1. `data_preprocessing.R` reads the official RNA matrix and MULTI-seq table.
2. It scores an erythroid/K562 marker program and selects the two highest
   scoring MULTI-seq tag groups. The selected tags and scores are recorded in
   `SatpathyK562.multiseq_tag_summary.tsv` and the QC summary.
3. It translates the selected GEX barcodes to the paired ATAC barcodes using
   the line-aligned 10x `737K-arc-v1` GEX/ATAC whitelists, then counts the
   released GRCh38 peaks.  `fragment_qc.cpp` (or the Python fallback)
   calculates fragment count, nucleosome signal and TSS enrichment in one
   streaming pass; this avoids three full scans of the 7.8-GB compressed file.
4. It applies the K562 thresholds used in the supplied QC table:
   `nCount_RNA > 1000`, `nCount_ATAC > 200`, `NucleosomeRatio < 2`, and
   `TSSEnrichment > 2`.
5. `gene_peak_preprocessing.R` creates the +/-500-kb candidate
   `SatpathyK562.gene_peak_df.rds` used by EpiLink.

Run from the repository root:

```bash
bash data/SatpathyK562/download_official.sh
SATPATHY_PEAK_LIMIT=0 Rscript data/SatpathyK562/data_preprocessing.R
Rscript data/SatpathyK562/gene_peak_preprocessing.R
Rscript data/SatpathyK562/plot_celltype_umap.R
```

The generated files are `SatpathyK562.rna.rds`, `SatpathyK562.atac.rds`,
`SatpathyK562.meta.rds`, and `SatpathyK562.gene_peak_df.rds`. The current
benchmark run retains 6,071 paired QC-passing cells and all 189,805 released
peaks. On a smaller node, setting `SATPATHY_PEAK_LIMIT=30000` selects a
reproducible evenly spaced peak subset. The raw
fragments file is approximately 7.8 GB and is required to calculate ATAC
counts, nucleosome signal and TSS enrichment; it should not be committed to
git.  `SatpathyK562.qc_summary.tsv` records the thresholds and summary values;
`SatpathyK562.celltype_umap.pdf` is the RNA UMAP coloured by `newCT`.
Do not use `SATPATHY_SKIP_FRAGMENT_QC=1` for benchmark results; that option is
only a debugging mode that leaves nucleosome/TSS metrics missing.
