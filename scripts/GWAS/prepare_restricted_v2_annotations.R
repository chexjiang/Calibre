#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
  library(rtracklayer)
})

args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[[1]] else ".", mustWork = TRUE)
gwas <- file.path(root, "results", "validation", "GWAS")
version <- file.path(gwas, "SCENT_style", "testing", "restricted_v2")
methods <- fread(file.path(version, "methods.tsv"))
datasets <- fread(file.path(gwas, "SCENT_style", "config", "datasets.tsv"))

# Requested cell-type restriction. Other datasets retain their existing
# inventory, so NEAT/Brain/K562 remain available in the same version.
keep <- list(
  pbmc = c("T", "Bcell", "Myeloid"),
  bmmc = c("T", "Bcell", "Myeloid"),
  Dogma_CTRL = "T",
  Dogma_STIM = "T",
  SHARE = "LCL",
  Arthritis = c("T", "B", "Myeloid")
)
sel <- rep(TRUE, nrow(datasets))
for (nm in names(keep)) {
  sel[datasets$dataset == nm] <- datasets$celltype[datasets$dataset == nm] %chin% keep[[nm]]
}
datasets <- datasets[sel]

chain_19_38 <- file.path(gwas, "reference", "liftover", "hg19ToHg38.over.chain.gz")
chain_38_19 <- file.path(gwas, "sLDSC", "reference", "raw", "hg38ToHg19.over.chain.gz")
read_chain <- function(path) {
  tmp <- tempfile(fileext=".chain")
  writeLines(readLines(gzfile(path), warn=FALSE), tmp)
  import.chain(tmp)
}
chain19to38 <- read_chain(chain_19_38)
chain38to19 <- read_chain(chain_38_19)

safe_id <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
parse_peaks <- function(x) {
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) return(data.table(chr=character(), start=integer(), end=integer()))
  z <- tstrsplit(unique(x), "-", fixed=TRUE)
  if (length(z) != 3L) stop("Unexpected peak format")
  ans <- data.table(chr=z[[1]], start=as.integer(z[[2]]), end=as.integer(z[[3]]))
  ans[grepl("^chr([1-9]|1[0-9]|2[0-2]|X)$", chr) & !is.na(start) & !is.na(end) & start <= end]
}
reduce_bed <- function(x) {
  if (!nrow(x)) return(x)
  gr <- reduce(GRanges(x$chr, IRanges(x$start, x$end)), ignore.strand=TRUE)
  data.table(chr=as.character(seqnames(gr)), start=start(gr), end=end(gr))
}
lift_bed <- function(x, chain) {
  if (!nrow(x)) return(x)
  y <- unlist(liftOver(GRanges(x$chr, IRanges(x$start, x$end)), chain), use.names=FALSE)
  if (!length(y)) return(data.table(chr=character(), start=integer(), end=integer()))
  reduce_bed(data.table(chr=as.character(seqnames(y)), start=start(y), end=end(y)))
}
write_bed <- function(x, path) {
  dir.create(dirname(path), recursive=TRUE, showWarnings=FALSE)
  y <- copy(x); y[, start0 := pmax(start - 1L, 0L)]
  fwrite(y[, .(chr, start0, end)], path, sep="\t", col.names=FALSE)
}
read_bed <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0L)
    return(data.table(chr=character(), start=integer(), end=integer()))
  fread(path, header=FALSE, col.names=c("chr", "start0", "end"))[, .(
    chr, start=start0 + 1L, end
  )]
}

ct_scent <- file.path(version, "celltype", "SCENT", "annotations")
ct_sldsc <- file.path(gwas, "sLDSC", "testing", "restricted_v2", "celltype", "sLDSC", "annotations")
pooled_scent <- file.path(version, "pooled", "SCENT", "annotations")
pooled_sldsc <- file.path(gwas, "sLDSC", "testing", "restricted_v2", "pooled", "sLDSC", "annotations")
dir.create(ct_scent, recursive=TRUE, showWarnings=FALSE)
dir.create(ct_sldsc, recursive=TRUE, showWarnings=FALSE)
dir.create(pooled_scent, recursive=TRUE, showWarnings=FALSE)
dir.create(pooled_sldsc, recursive=TRUE, showWarnings=FALSE)

ct_rows <- list(); k <- 0L
for (i in seq_len(nrow(datasets))) for (j in seq_len(nrow(methods))) {
  d <- datasets[i]; m <- methods[j]
  input <- file.path(root, "results", "real", d$dataset, d$celltype,
                     m$source_method, m$source_file)
  if (!file.exists(input)) next
  dt <- fread(input, select="peak", showProgress=FALSE)
  native <- reduce_bed(parse_peaks(dt$peak))
  hg38 <- if (d$build == "hg38") native else lift_bed(native, chain19to38)
  hg19 <- if (d$build == "hg19") native else lift_bed(native, chain38to19)
  id <- paste(safe_id(d$dataset), safe_id(d$celltype), safe_id(m$method), sep="__")
  write_bed(hg38, file.path(ct_scent, paste0(id, ".hg38.bed")))
  write_bed(hg19, file.path(ct_sldsc, paste0(id, ".hg19.bed")))
  k <- k + 1L
  ct_rows[[k]] <- data.table(
    annotation_id=id, dataset=d$dataset, celltype=d$celltype,
    panel=d$panel, interpretation=d$interpretation, method=m$method,
    method_label=m$label, source_method=m$source_method,
    source_file=m$source_file, source_build=d$build,
    n_pairs=nrow(dt), n_native_intervals=nrow(native),
    n_hg38_intervals=nrow(hg38), n_hg19_intervals=nrow(hg19),
    lift_fraction_hg38=if (nrow(native)) nrow(hg38)/nrow(native) else NA_real_,
    lift_fraction_hg19=if (nrow(native)) nrow(hg19)/nrow(native) else NA_real_
  )
}
ct_manifest <- rbindlist(ct_rows, fill=TRUE)
fwrite(ct_manifest, file.path(version, "celltype", "annotation_manifest.tsv"), sep="\t")
fwrite(ct_manifest, file.path(gwas, "sLDSC", "testing", "restricted_v2", "celltype", "annotation_manifest.tsv"), sep="\t")

pooled_rows <- list(); k <- 0L
groups <- unique(ct_manifest[, .(dataset, panel)])
for (gi in seq_len(nrow(groups))) for (m in unique(methods$method)) {
  g <- groups[gi]
  z <- ct_manifest[dataset == g$dataset & panel == g$panel & method == m]
  if (!nrow(z)) next
  hg38 <- rbindlist(lapply(z$annotation_id, function(id)
    read_bed(file.path(ct_scent, paste0(id, ".hg38.bed")))), fill=TRUE)
  hg19 <- rbindlist(lapply(z$annotation_id, function(id)
    read_bed(file.path(ct_sldsc, paste0(id, ".hg19.bed")))), fill=TRUE)
  hg38 <- reduce_bed(hg38); hg19 <- reduce_bed(hg19)
  id <- paste(safe_id(g$dataset), "ALL", safe_id(g$panel), safe_id(m), sep="__")
  write_bed(hg38, file.path(pooled_scent, paste0(id, ".hg38.bed")))
  write_bed(hg19, file.path(pooled_sldsc, paste0(id, ".hg19.bed")))
  k <- k + 1L
  pooled_rows[[k]] <- data.table(
    annotation_id=id, dataset=g$dataset, celltype="ALL", panel=g$panel,
    interpretation=if (all(z$interpretation == "primary")) "primary" else "mixed",
    method=m, method_label=methods[method == m, label][1],
    source_method=methods[method == m, source_method][1],
    source_build="mixed", n_source_celltypes=uniqueN(z$celltype),
    n_source_annotations=nrow(z), n_hg38_intervals=nrow(hg38), n_hg19_intervals=nrow(hg19)
  )
}
pooled_manifest <- rbindlist(pooled_rows, fill=TRUE)
fwrite(pooled_manifest, file.path(version, "pooled", "annotation_manifest.tsv"), sep="\t")
fwrite(pooled_manifest, file.path(gwas, "sLDSC", "testing", "restricted_v2", "pooled", "annotation_manifest.tsv"), sep="\t")
message("Cell-type annotations: ", nrow(ct_manifest), "; pooled annotations: ", nrow(pooled_manifest))
