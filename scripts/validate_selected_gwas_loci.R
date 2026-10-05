#!/usr/bin/env Rscript

# Reproduce exact interval overlaps between selected EpiLink GWAS peaks and
# locally available CRISPR enhancer-gene benchmarks (GRCh38).

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[[1L]] else ".", mustWork = TRUE)
out_dir <- file.path(root, "results", "validation", "functional", "GWAS_loci")

selected <- fread(file.path(out_dir, "selected_loci_evidence.tsv"))
crispr <- fread(file.path(root, "reference", "CRISPR", "crispr_hg38_positive.tsv"))

peak_bits <- tstrsplit(selected$peak, "-", fixed = TRUE)
query <- selected[, .(locus_id, variant, gene, peak)]
query[, `:=`(chr = peak_bits[[1L]], peak_start = as.integer(peak_bits[[2L]]),
             peak_end = as.integer(peak_bits[[3L]]))]
query[, query_id := .I]

setnames(crispr, c("start", "end"), c("crispr_start", "crispr_end"))
query_iv <- query[, .(chr, start = peak_start, end = peak_end, query_id)]
crispr_iv <- crispr[, .(chr, start = crispr_start, end = crispr_end,
                        crispr_gene = gene, regulated)]
setkey(query_iv, chr, start, end)
setkey(crispr_iv, chr, start, end)

hits <- foverlaps(query_iv, crispr_iv, type = "any", nomatch = 0L)
if (nrow(hits)) {
  hits <- merge(hits, query, by = "query_id", suffixes = c("_crispr", "_epilink"))
  if ("chr_crispr" %chin% names(hits)) setnames(hits, "chr_crispr", "chr")
  if ("chr.x" %chin% names(hits)) setnames(hits, "chr.x", "chr")
  hits[, same_gene := gene == crispr_gene]
  hits[, overlap_bp := pmax(0L, pmin(end, peak_end) - pmax(start, peak_start) + 1L)]
  hits <- hits[, .(locus_id, variant, gene, peak, chr, peak_start, peak_end,
                   crispr_gene, crispr_start = start, crispr_end = end,
                   regulated, same_gene, overlap_bp)]
}

fwrite(hits, file.path(out_dir, "selected_loci_crispr_overlaps.tsv"), sep = "\t")
message("Wrote ", nrow(hits), " CRISPR interval overlaps; ",
        hits[, sum(same_gene & regulated)], " regulate the nominated target gene")
