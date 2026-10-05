Arthritis tissue benchmark preparation

Source: GEO GSE243917. The complete raw gene/peak H5AD matrices are already post-QC.
The paper thresholds are recorded in `Arthritis.qc_thresholds.tsv`.

The complete GEO object has 654 mural cells; these are excluded to reproduce the
published 30,893-cell SCENT arthritis panel. Bcell+plasma and Tcell+NK are merged.
`Arthritis.atac.rds` is a peak-by-cell count matrix; SCENT should binarize it
method-side. Candidate pairs use hg38 gene bodies +/-500 kb and the shared peak set.
