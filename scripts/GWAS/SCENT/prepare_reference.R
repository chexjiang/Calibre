#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
  library(rtracklayer)
})

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)
gwas <- file.path(root, "results", "validation", "GWAS")
finngen_ref <- file.path(gwas, "reference", "FinnGen")
common_ref <- file.path(gwas, "reference", "1000G")
chain_ref <- file.path(gwas, "reference", "liftover")
ldref <- file.path(gwas, "sLDSC", "reference", "processed")
dir.create(finngen_ref, recursive = TRUE, showWarnings = FALSE)
dir.create(common_ref, recursive = TRUE, showWarnings = FALSE)

traits <- unique(fread(file.path(gwas, "SCENT_style", "config", "traits.tsv"))[, .(trait, label, phenotype, cases, controls)])
fm <- list()
for (tr in traits$trait) {
  f <- file.path(finngen_ref, paste0("finngen_R12_", tr, ".SUSIE.snp.filter.tsv"))
  x <- fread(f)
  setnames(x, "cs_specific_prob", "pip")
  x[, chromosome := sub("^chr", "", chromosome)]
  x <- x[, .(pip = max(pip, na.rm = TRUE), region = region[which.max(pip)], maf = maf[which.max(pip)]),
         by = .(trait, chromosome, position, allele1, allele2)]
  x[, chr := paste0("chr", chromosome)]
  fm[[tr]] <- x
}
finemap <- rbindlist(fm, use.names = TRUE, fill = TRUE)
saveRDS(finemap, file.path(finngen_ref, "finngen_R12_susie_all.rds"), compress = "xz")
fwrite(finemap[pip >= 0.2], file.path(finngen_ref, "finngen_R12_PIPge0.2.tsv.gz"), sep = "\t")

find_one <- function(pattern) {
  hits <- list.files(ldref, pattern = pattern, recursive = TRUE, full.names = TRUE)
  if (!length(hits)) stop("Missing extracted LDSC reference matching: ", pattern)
  hits[[1]]
}

import_chain_gz <- function(path) {
  tmp <- tempfile(fileext = ".over.chain")
  in_con <- gzfile(path, open = "rt")
  lines <- readLines(in_con, warn = FALSE)
  close(in_con)
  writeLines(lines, tmp)
  import.chain(tmp)
}
chain <- import_chain_gz(file.path(chain_ref, "hg19ToHg38.over.chain.gz"))
common <- vector("list", 22L)
for (chr in 1:22) {
  bim <- fread(find_one(paste0("1000G\\.(EUR\\.QC|mac5eur)\\.", chr, "\\.bim$")), header = FALSE,
               col.names = c("CHR", "SNP", "CM", "BP", "A1", "A2"))
  frq <- fread(find_one(paste0("1000G\\.EUR\\.QC\\.", chr, "\\.frq$")))
  keep <- frq[MAF > 0.01, .(SNP, MAF)]
  z <- bim[keep, on = "SNP", nomatch = 0L]
  gr <- GRanges(paste0("chr", chr), IRanges(z$BP, z$BP), SNP = z$SNP, MAF = z$MAF)
  lifted <- liftOver(gr, chain)
  ok <- lengths(lifted) == 1L
  one <- unlist(lifted[ok], use.names = FALSE)
  common[[chr]] <- data.table(
    chr = as.character(seqnames(one)), pos = start(one),
    SNP = mcols(gr)$SNP[ok], MAF = mcols(gr)$MAF[ok]
  )
  message("chr", chr, ": ", nrow(common[[chr]]), " common EUR variants lifted to hg38")
}
common <- unique(rbindlist(common), by = c("chr", "pos"))
saveRDS(common, file.path(common_ref, "1000G_EUR_MAFgt0.01.hg38.rds"), compress = FALSE)

manifest <- data.table(
  resource = c("FinnGen R12 SuSiE", "1000G Phase3 EUR common variants"),
  build = c("hg38", "hg38 (lifted from hg19)"),
  threshold = c("retained all; analysis uses PIP >= 0.2", "MAF > 0.01"),
  n = c(nrow(finemap), nrow(common))
)
fwrite(manifest, file.path(finngen_ref, "reference_manifest.tsv"), sep = "\t")
