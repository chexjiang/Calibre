#!/usr/bin/env Rscript

file_argument <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))])
project_root <- normalizePath(file.path(dirname(file_argument), ".."), mustWork = TRUE)
setwd(project_root)

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("data.table is required", call. = FALSE)
}
if (!nzchar(Sys.which("bedtools"))) {
  stop("bedtools is required to prepare CRISPR_new", call. = FALSE)
}

reference_dir <- file.path(project_root, "reference", "CRISPR_new")
element_file <- file.path(reference_dir, "EPCrisprBenchmark_ensemble_data_GRCh38.tsv.gz")
snp_file <- file.path(reference_dir, "SNPs.bed.gz")
output_file <- file.path(reference_dir, "crispr_new_snp_gene_benchmark.tsv.gz")
summary_file <- file.path(reference_dir, "benchmark_summary.tsv")

if (!file.exists(element_file)) stop("Missing ", element_file, call. = FALSE)
if (!file.exists(snp_file)) stop("Missing ", snp_file, call. = FALSE)

elements <- data.table::fread(element_file, showProgress = FALSE)
required <- c("chrom", "chromStart", "chromEnd", "measuredGeneSymbol",
              "startTSS", "Regulated", "ValidConnection", "CellType", "Dataset")
if (!all(required %in% names(elements))) {
  stop("CRISPR benchmark lacks: ",
       paste(setdiff(required, names(elements)), collapse = ", "), call. = FALSE)
}
# A missing Regulated value has no positive/negative experimental label and
# therefore cannot be used in a binary benchmark.
elements <- elements[
  as.character(ValidConnection) == "TRUE" & CellType == "K562" &
    !is.na(Regulated) & !is.na(startTSS)
]
elements[, element_id := .I]

element_bed <- tempfile("crispr_elements_", fileext = ".bed")
overlap_file <- tempfile("crispr_snp_overlap_", fileext = ".tsv")
bedtools_log <- tempfile("crispr_bedtools_", fileext = ".log")
on.exit(unlink(c(element_bed, overlap_file, bedtools_log)), add = TRUE)
data.table::fwrite(
  elements[, .(chrom, chromStart, chromEnd, element_id,
               measuredGeneSymbol, Regulated, startTSS, Dataset)],
  element_bed, sep = "\t", col.names = FALSE
)

status <- system2(
  Sys.which("bedtools"),
  c("intersect", "-a", element_bed, "-b", snp_file, "-wa", "-wb"),
  stdout = overlap_file, stderr = bedtools_log
)
if (!identical(status, 0L)) {
  details <- if (file.exists(bedtools_log)) {
    paste(readLines(bedtools_log, warn = FALSE), collapse = "\n")
  } else ""
  stop("bedtools intersect failed", if (nzchar(details)) paste0(":\n", details),
       call. = FALSE)
}

overlap <- data.table::fread(
  overlap_file, header = FALSE, showProgress = FALSE,
  col.names = c("element_chr", "element_start", "element_end", "element_id",
                "gene", "regulated", "tss", "dataset",
                "snp_chr", "snp_start", "snp_end")
)
overlap[, positive_row := toupper(as.character(regulated)) == "TRUE"]
overlap[, snp := paste(snp_chr, snp_start, snp_end, sep = "-")]
overlap[, distance_to_tss := abs((snp_start + 1L) - as.integer(tss))]

# A SNP-gene pair can arise from overlapping tested elements. A validated
# regulatory result takes precedence over a non-significant result.
benchmark <- overlap[, .(
  chr = snp_chr[[1L]],
  snp_start = min(snp_start),
  snp_end = max(snp_end),
  positive = any(positive_row),
  distance_to_tss = min(distance_to_tss),
  n_tested_elements = data.table::uniqueN(element_id),
  datasets = paste(sort(unique(dataset)), collapse = ";")
), by = .(snp, gene)]
data.table::setorder(benchmark, chr, snp_start, gene)
data.table::fwrite(benchmark, output_file, sep = "\t", compress = "gzip")

summary <- data.frame(
  stage = c("K562 element-gene tests", "SNP-gene benchmark", "positive", "negative"),
  n = c(nrow(elements), nrow(benchmark), sum(benchmark$positive),
        sum(!benchmark$positive))
)
utils::write.table(summary, summary_file, sep = "\t", quote = FALSE,
                   row.names = FALSE)
message("Wrote ", nrow(benchmark), " SNP-gene benchmark rows (",
        sum(benchmark$positive), " positive; ", sum(!benchmark$positive),
        " negative) to ", output_file)
