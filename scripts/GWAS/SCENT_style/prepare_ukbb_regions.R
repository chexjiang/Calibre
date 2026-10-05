#!/usr/bin/env Rscript

# Build the UKBB hg38 locus table separately from the causal-variant table.
# This keeps memory bounded and uses the complete regions BED shipped with the
# release, rather than only regions represented by PIP >= 0.2 variants.

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
  library(rtracklayer)
})

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)
gwas <- file.path(root, "results", "validation", "GWAS")
raw <- file.path(gwas, "reference", "UKBB")
out_dir <- file.path(gwas, "reference", "UKBB")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

import_chain_gz <- function(path) {
  tmp <- tempfile(fileext = ".over.chain")
  con <- gzfile(path, open = "rt")
  lines <- readLines(con, warn = FALSE)
  close(con)
  writeLines(lines, tmp)
  import.chain(tmp)
}
chain <- import_chain_gz(file.path(gwas, "reference", "liftover", "hg19ToHg38.over.chain.gz"))

regions_bed <- file.path(raw, "UKBB_94traits_release1_regions.bed.gz")
if (!file.exists(regions_bed)) stop("Missing UKBB regions BED: ", regions_bed)
region_cmd <- paste("gzip -dc", shQuote(regions_bed),
                    "| awk 'BEGIN{FS=\"\\t\"; OFS=\"\\t\"} {k=$5 OFS $6; if (!seen[k]++) print $5, $6}'")
regions <- fread(cmd = region_cmd, header = FALSE, col.names = c("trait", "region"), showProgress = TRUE)
regions <- regions[!is.na(region) & nzchar(region)]
z <- tstrsplit(regions$region, "[:-]")
gr19 <- GRanges(z[[1]], IRanges(as.integer(z[[2]]), as.integer(z[[3]])))
lifted <- liftOver(gr19, chain)
ok <- lengths(lifted) > 0L
idx <- rep(which(ok), lengths(lifted[ok]))
gr38 <- unlist(lifted[ok], use.names = FALSE)
out <- data.table(
  trait = regions$trait[idx],
  chr = as.character(seqnames(gr38)),
  start = start(gr38), end = end(gr38)
)

reduce_one <- function(x) {
  g <- reduce(GRanges(x$chr, IRanges(x$start, x$end)), ignore.strand = TRUE)
  data.table(trait = x$trait[1], chr = as.character(seqnames(g)), start = start(g), end = end(g))
}
out <- rbindlist(lapply(split(out, out$trait), reduce_one), fill = TRUE)
out <- out[!(chr == "chr6" & end >= 25000000L & start <= 36000000L)]
fwrite(out, file.path(out_dir, "ukbb_94traits_regions.hg38.tsv.gz"), sep = "\t")
message("Wrote ", nrow(out), " reduced UKBB hg38 trait-locus intervals for ", uniqueN(out$trait), " traits")
