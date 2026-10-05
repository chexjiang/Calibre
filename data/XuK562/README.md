# Xu K562 10x Multiome

This directory contains the paired RNA/ATAC matrices and metadata for the
Xu K562 benchmark dataset. The original study is ArrayExpress accession
`E-MTAB-11264` (K562 10x Multiome; ENA runs `ERR9847049` for ATAC and
`ERR9847050` for RNA). The matrices currently stored here were converted from
the public Cell2Net `K562.h5mu` object, which was generated with a Cell Ranger
ARC workflow. The source object has 6,508 paired cells; the stored matrices
contain the 750 cells passing the Xu K562 QC thresholds below.

Files:

- `XuK562.rna.rds`: gene-by-cell raw RNA counts, with gene symbols.
- `XuK562.atac.rds`: peak-by-cell raw ATAC counts, GRCh38 peaks in
  `chr-start-end` format.
- `XuK562.meta.rds`: paired-cell metadata (`barcodes`, `newCT = K562`,
  `nCount_RNA`, `nCount_ATAC`, `NucleosomeRatio`, `TSSEnrichment`, and
  `qc_pass`).
- `XuK562.gene_peak_df.rds`: +/-500-kb candidate gene-peak pairs generated
  from the project hg38 annotation.

QC applied before saving the matrices: `nCount_RNA > 1000`,
`nCount_ATAC > 200`, `NucleosomeRatio < 2`, and `TSSEnrichment > 2`.

The Dorans and Price manuscript reports 7,861 cells and 132,736 peaks after
its own Seurat/ArchR QC and MACS2 peak-calling. The public h5mu object used
here contains 6,508 cells and 135,964 peaks before QC, so these files are
format-compatible Xu K562 inputs but are not a byte-for-byte recreation of
the authors' intermediate matrices.

To recreate the conversion from the downloaded h5mu object:

```bash
curl -L 'https://www.dropbox.com/scl/fi/8kkyyinosak213ksmp7p9/K562.h5mu?dl=1&rlkey=4nlcc6zqf2y5lv2at6cw49ygg&st=5hs15lgp' -o /tmp/K562.h5mu
Rscript data/XuK562/import_h5mu.R
Rscript data/XuK562/data_preprocessing.R
```

`import_h5mu.R` expects `/tmp/K562.h5mu` by default; override it with
`XUK562_H5MU=/path/to/K562.h5mu`. The h5mu download link and processing
description are documented in the [Cell2Net K562 tutorial](https://pinellolab.github.io/cell2net/docs/latest/tutorials/k562/01_pretrain_seq_encoder.html);
the authoritative raw-data source remains [ArrayExpress E-MTAB-11264](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-11264).
