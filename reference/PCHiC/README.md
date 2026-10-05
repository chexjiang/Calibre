# PCHiC reference

This directory contains the Javierre et al. promoter-capture Hi-C reference
used by the real-data validation workflow. Coordinates in the peak matrix are
lifted from hg19 to hg38 with `hg19ToHg38.over.chain`.

The active implementation is `scripts/PCHiC_enrichment.R`. It uses the other-end
fragment, exact `baitName` matching, CHiCAGO score `> 5` (the executable
scMultiMap procedure), and a one-sided Fisher exact test. The broad EpiLink
groups are mapped to `nB/tB`, `Mon`, and the available CD4/CD8 T-cell columns;
there is no NK-specific PCHiC column, so `Tnk` is an explicitly documented
proxy comparison.

For SHARE-seq LCL, the runner uses `nB/tB` as a primary B-cell proxy and lifts
the hg19 SHARE candidate peak ranges to hg38 at validation time. The original
SHARE `gene--peak` keys are retained.
