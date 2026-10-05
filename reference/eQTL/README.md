# eQTL validation references

This directory contains the immutable raw inputs and normalized files used by
the two EpiLink eQTL validation pipelines.

## Primary fine-mapping reference

The primary reference is release-7 eQTL Catalogue SuSiE fine-mapping for the
OneK1K study (`QTS000038`).  Eight gene-expression (`ge`) datasets are grouped
into B cell, myeloid, and T/NK benchmarks in `config/eqtl_validation.yml`.
Within a broad group, duplicate variant-gene pairs are collapsed using their
maximum PIP.  The source files and metadata are public:

- https://www.ebi.ac.uk/eqtl/Data_access/
- https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/susie/QTS000038/

The processed table uses GRCh38, one-based variant coordinates, HGNC symbols,
and retains all credible-set members so that either PIP >= 0.20 (SCENT-style)
or PIP >= 0.50 (ranking-style) can be applied without reprocessing.

## Background variants

`raw/pgboost/1000G_MAC5_SNPs.hg38.bed.gz` is the 1000 Genomes SNP set with
minor allele count >5 distributed in the pgBoost repository.  BED coordinates
are converted from zero-based half-open to one-based closed coordinates during
preparation.  The accompanying GTEx training labels are archived for
provenance but are not used to train or validate EpiLink.

- https://github.com/elizabethdorans/pgBoost

## Reproducibility

Run `scripts/eQTL_shared/download_references.sh` to re-download/copy the exact
inputs, then `Rscript scripts/eQTL_shared/prepare_reference.R`.  Checksums and
row counts are written to `shared/processed/reference_manifest.tsv`.

