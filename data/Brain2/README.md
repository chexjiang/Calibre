Brain2 benchmark preparation

Source: GEO GSE214979, Anderson et al. (2023), postmortem DLPFC snRNA+snATAC multiome.
The filtered 10x H5 contains 105,332 nuclei, 36,601 genes and 150,614 peaks.
The supplied cell metadata contains RNA/ATAC depth, TSS and nucleosome QC metrics.
The source has 8 annotated cell types; this benchmark export retains cell types with >500 cells, leaving 6 types.
The source-level 8-type counts are in Brain2.source_celltype_counts.tsv.
Candidate gene_peak_df uses hg38 gene bodies +/-500 kb and the filtered peak set.
The scMultiMap paper reports 380,154 peaks in its summary table; this differs from the 150,614 peak features in the GEO filtered H5.
The 59 GB GEO fragment file was not downloaded; the supplied metadata is used for QC and the filtered count matrix for benchmarking.
