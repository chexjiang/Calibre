# H3K27ac HiChIP reference

`41588_2017_BFng3963_MOESM4_ESM.xlsx` contains the hg19 loop coordinates from
the H3K27ac HiChIP study used by SCENT's orthogonal consistency analysis.
`hg19ToHg38.over.chain` is used by `scripts/run_hichip_validation.R`.

The implementation follows the paper's enhancer--gene construction: each loop
anchor is lifted to hg38, gene promoters are defined as 1 kb upstream and 1 kb
downstream of hg38 gene TSSs, and a loop is retained when one anchor overlaps a
promoter and the other anchor is the enhancer interval. Naive, Th17, and Treg
links are pooled as the Tnk reference. GM12878 is used as the Bcell reference;
it is a B-lymphoblastoid cell-line proxy rather than a primary B-cell map.

For SHARE-seq LCL, `GM12878 H3K27ac Loops` is used and the hg19 SHARE candidate
peak ranges are lifted to hg38 during validation while preserving the original
pair keys.
