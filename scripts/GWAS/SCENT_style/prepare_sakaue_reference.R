#!/usr/bin/env Rscript

# Prepare the reference sets used by the Sakaue et al.-style SCENT analysis.
#
# The analysis uses:
#   causal variants: fine-mapped variants with PIP >= 0.2
#   denominator: 1000 Genomes EUR variants with MAF > 1% in the fine-mapped loci
#   score: (causal in annotation / common in annotation) /
#          (causal in loci / common in loci)
#
# FinnGen is already in hg38 in this project. UK Biobank release 1.1 is in
# hg19 and is lifted to hg38 here so both references use the same annotation BEDs.

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
  library(rtracklayer)
})

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)
gwas <- file.path(root, "results", "validation", "GWAS")
base <- file.path(gwas, "SCENT_style", "UKBB", "sakaue_style")
ukbb_ref <- file.path(gwas, "reference", "UKBB")
finngen_ref <- file.path(gwas, "reference", "FinnGen")
common_ref <- file.path(gwas, "reference", "1000G")
chain_ref <- file.path(gwas, "reference", "liftover")
dir.create(ukbb_ref, recursive = TRUE, showWarnings = FALSE)
dir.create(finngen_ref, recursive = TRUE, showWarnings = FALSE)

common <- readRDS(file.path(common_ref, "1000G_EUR_MAFgt0.01.hg38.rds"))
common <- unique(common[, .(chr = as.character(chr), pos = as.integer(pos), SNP, MAF)], by = c("chr", "pos"))

import_chain_gz <- function(path) {
  tmp <- tempfile(fileext = ".over.chain")
  con <- gzfile(path, open = "rt")
  lines <- readLines(con, warn = FALSE)
  close(con)
  writeLines(lines, tmp)
  import.chain(tmp)
}

chain_19_38 <- import_chain_gz(file.path(chain_ref, "hg19ToHg38.over.chain.gz"))

empty_intervals <- function() data.table(chr = character(), start = integer(), end = integer())
empty_points <- function() data.table(chr = character(), pos = integer())

drop_mhc <- function(x) {
  if (!nrow(x)) return(x)
  x[!(chr == "chr6" & pos >= 25000000L & pos <= 36000000L)]
}

parse_regions <- function(x) {
  x <- unique(x[!is.na(x) & nzchar(x)])
  if (!length(x)) return(empty_intervals())
  z <- tstrsplit(x, "[:-]")
  reduce(GRanges(z[[1]], IRanges(as.integer(z[[2]]), as.integer(z[[3]]))), ignore.strand = TRUE) |>
    (
      function(gr) data.table(chr = as.character(seqnames(gr)), start = start(gr), end = end(gr))
    )()
}

lift_points <- function(x) {
  if (!nrow(x)) return(empty_points())
  gr <- GRanges(x$chr, IRanges(x$pos, x$pos))
  lifted <- liftOver(gr, chain_19_38)
  ok <- lengths(lifted) == 1L
  if (!any(ok)) return(empty_points())
  one <- unlist(lifted[ok], use.names = FALSE)
  unique(data.table(chr = as.character(seqnames(one)), pos = as.integer(start(one))))
}

lift_intervals <- function(x) {
  if (!nrow(x)) return(empty_intervals())
  gr <- GRanges(x$chr, IRanges(x$start, x$end))
  lifted <- unlist(liftOver(gr, chain_19_38), use.names = FALSE)
  if (!length(lifted)) return(empty_intervals())
  out <- GRanges(seqnames(lifted), ranges(lifted))
  out <- reduce(out, ignore.strand = TRUE)
  data.table(chr = as.character(seqnames(out)), start = start(out), end = end(out))
}

overlap_mask <- function(points, intervals) {
  hit <- rep(FALSE, nrow(points))
  if (!nrow(points) || !nrow(intervals)) return(hit)
  by_chr <- split(seq_len(nrow(points)), points$chr)
  for (cc in intersect(names(by_chr), unique(intervals$chr))) {
    rows <- by_chr[[cc]]
    iv <- intervals[chr == cc][order(start)]
    idx <- findInterval(points$pos[rows], iv$start)
    ok <- idx > 0L
    local <- rep(FALSE, length(rows))
    local[ok] <- points$pos[rows][ok] <= iv$end[idx[ok]]
    hit[rows] <- local
  }
  hit
}

message("Preparing FinnGen reference")
fg <- readRDS(file.path(finngen_ref, "finngen_R12_susie_all.rds"))
fg_causal <- fg[pip >= 0.2, .(trait, chr = as.character(chr), pos = as.integer(position), pip)]
fg_causal <- drop_mhc(unique(fg_causal, by = c("trait", "chr", "pos")))
fg_regions <- fg[, .(regions = list(parse_regions(region))), by = trait]
fg_regions[, regions := lapply(regions, function(x) {
  if (!nrow(x)) return(x)
  x[!(chr == "chr6" & end >= 25000000L & start <= 36000000L)]
})]
fg_trait_names <- sort(unique(fg$trait))

fg_ref <- list(
  source = "FinnGen R12 project subset",
  build = "hg38",
  pip_threshold = 0.2,
  causal = fg_causal,
  regions = setNames(fg_regions$regions, fg_regions$trait),
  trait_names = fg_trait_names,
  common = common
)
saveRDS(fg_ref, file.path(finngen_ref, "finngen_sakaue_reference.rds"), compress = "xz")

# Release the large common-variant table before reading the UKBB stream.  It
# is reloaded only after UKBB causal/locus coordinates have been lifted.
rm(fg, fg_ref, fg_causal, fg_regions, common)
gc(verbose = FALSE)

message("Preparing UK Biobank reference")
ukbb_dir <- file.path(ukbb_ref)
ukbb_bed <- file.path(ukbb_dir, "UKBB_94traits_release1.bed.gz")
if (!file.exists(ukbb_bed)) stop("Missing UKBB primary BED: ", ukbb_bed)
uk_cols <- c(
  "chromosome", "start0", "end0", "variant", "rsid", "allele1", "allele2",
  "minorallele", "cohort", "model_marginal", "method", "trait", "region",
  "maf", "beta_marginal", "se_marginal", "chisq_marginal", "pip", "cs_id",
  "beta_posterior", "sd_posterior", "LD_HWE", "LD_SV"
)
# The primary BED is large.  Stream-filter on PIP before fread so the R
# process never materializes the full 94-trait file in memory.
uk_cmd <- paste("gzip -dc", shQuote(ukbb_bed), "| awk '$18 >= 0.2'")
uk <- fread(cmd = uk_cmd, header = FALSE, col.names = uk_cols, showProgress = TRUE)
uk <- uk[!is.na(pip) & pip >= 0.2]
uk[, `:=`(chr = as.character(chromosome), pos19 = as.integer(start0) + 1L)]
uk <- uk[!(chr == "chr6" & pos19 >= 25000000L & pos19 <= 36000000L)]

uk_trait_names <- sort(unique(uk$trait))
uk_causal <- vector("list", length(uk_trait_names))
uk_regions <- vector("list", length(uk_trait_names))
names(uk_causal) <- names(uk_regions) <- uk_trait_names
region_dir <- file.path(ukbb_ref, "ukbb_regions_hg38")
dir.create(region_dir, recursive = TRUE, showWarnings = FALSE)
for (tr in uk_trait_names) {
  x <- uk[trait == tr]
  points19 <- unique(x[, .(chr, pos = pos19)])
  causal38 <- lift_points(points19)
  causal38[, trait := tr]
  uk_causal[[tr]] <- causal38
  reg19 <- parse_regions(x$region)
  reg19 <- reg19[!(chr == "chr6" & end >= 25000000L & start <= 36000000L)]
  regions38 <- lift_intervals(reg19)
  uk_regions[[tr]] <- regions38
  fwrite(cbind(regions38, trait = tr),
         file.path(region_dir, paste0(gsub("[^A-Za-z0-9]+", "_", tr), ".tsv.gz")),
         sep = "\t")
  message(tr, ": causal hg38 variants=", nrow(causal38), "; regions hg38=", nrow(regions38))
}
uk_causal_dt <- rbindlist(uk_causal, fill = TRUE)

# Keep UKBB as compact tabular files.  Regions are written one trait at a time
# to avoid materializing a large nested object during serialization.
fwrite(uk_causal_dt, file.path(ukbb_ref, "ukbb_94traits_causal.hg38.tsv.gz"), sep = "\t")

manifest <- rbindlist(list(
  data.table(source = "FinnGen", n_traits = length(fg_trait_names), n_causal = nrow(fg_causal), build = "hg38", pip_threshold = 0.2),
  data.table(source = "UKBB", n_traits = length(uk_trait_names), n_causal = nrow(uk_causal_dt), build = "hg38 lifted from hg19", pip_threshold = 0.2)
))
fwrite(manifest, file.path(ukbb_ref, "sakaue_reference_manifest.tsv"), sep = "\t")
message("Wrote Sakaue-style reference objects to ", ukbb_ref, " and ", finngen_ref)
