#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args)) args[[1L]] else "config/gtex_eqtl_validation.yml"
config <- read_yaml(config_path)
ref <- config$reference
dir.create(ref$processed_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(ref$caviar_file)) stop("Missing extracted CAVIAR file: ", ref$caviar_file)
if (!file.exists(ref$background)) stop("Missing 1000G EUR MAF>1% background: ", ref$background)

# Stream-filter the 286 MB gzip before fread materializes it. The paper uses
# PIP > 0.2 (strictly greater), and the >0.2 table also contains every record
# needed for the complementary PIP >0.5 ranked analysis.
cmd <- sprintf("gzip -dc %s | awk 'NR == 1 || $6 > %.17g'",
               shQuote(ref$caviar_file), as.numeric(ref$fine_map_pip_scent))
message("Reading GTEx CAVIAR records with PIP > ", ref$fine_map_pip_scent)
x <- fread(cmd = cmd, showProgress = TRUE)
setnames(x, c("TISSUE", "GENE", "eQTL", "CHROM", "POS", "Probability"),
         c("tissue", "ensembl_gene", "variant_id", "chr", "pos", "pip"))
x[, ensembl_gene := sub("\\..*$", "", ensembl_gene)]
x[, chr := ifelse(grepl("^chr", chr), chr, paste0("chr", chr))]
x[, `:=`(pos = as.integer(pos), pip = as.numeric(pip))]

gene_map <- as.data.table(readRDS(ref$ensembl_map))[, .(ensembl_gene = GENEID, gene = SYMBOL)]
x <- gene_map[x, on = "ensembl_gene", nomatch = 0L]
x <- x[!is.na(gene) & nzchar(gene) & !is.na(pos) & is.finite(pip)]
x <- unique(x[, .(tissue, gene, ensembl_gene, variant_id, chr,
                  start = pos, end = pos, pip)])

# SCENT only tests peaks whose center is within 500 kb of the target TSS. Keep
# the causal numerator/denominator in that same assessable cis universe.
tss <- as.data.table(readRDS(ref$tss))
x <- tss[x, on = .(gene, chr), nomatch = 0L]
x[, distance_to_tss := abs(start - tss)]
x <- x[distance_to_tss <= as.integer(ref$cis_window_bp)]

by_tissue <- x[, .(pip = max(pip), ensembl_gene = ensembl_gene[which.max(pip)],
                     distance_to_tss = distance_to_tss[which.max(pip)]),
                 by = .(tissue, gene, variant_id, chr, start, end)]
aggregated <- by_tissue[, .(
  pip = max(pip),
  n_tissues = uniqueN(tissue),
  tissues = paste(sort(unique(tissue)), collapse = ";"),
  distance_to_tss = distance_to_tss[which.max(pip)]
), by = .(gene, variant_id, chr, start, end)]

saveRDS(by_tissue, file.path(ref$processed_dir, "GTEx_v8_CAVIAR_PIPgt0.2_by_tissue.rds"),
        compress = FALSE)
saveRDS(aggregated, file.path(ref$processed_dir, "GTEx_v8_CAVIAR_PIPgt0.2_49tissues_aggregated.rds"),
        compress = FALSE)
fwrite(by_tissue, file.path(ref$processed_dir, "GTEx_v8_CAVIAR_PIPgt0.2_by_tissue.tsv.gz"), sep = "\t")
fwrite(aggregated, file.path(ref$processed_dir,
                            "GTEx_v8_CAVIAR_PIPgt0.2_49tissues_aggregated.tsv.gz"), sep = "\t")

tissues <- by_tissue[, .(
  n_variant_gene = .N,
  n_variants = uniqueN(variant_id),
  n_genes = uniqueN(gene),
  n_pip_gt_0.5 = sum(pip > as.numeric(ref$fine_map_pip_rank))
), by = tissue][order(tissue)]
fwrite(tissues, file.path(ref$processed_dir, "tissue_manifest.tsv"), sep = "\t")

background <- as.data.table(readRDS(ref$background))
audit <- data.table(
  item = c("CAVIAR tissue-specific variant-gene rows PIP>0.2 in +/-500kb",
           "CAVIAR aggregated unique variant-gene rows",
           "GTEx tissues", "Mapped genes", "1000G Phase3 EUR MAF>0.01 variants"),
  value = c(nrow(by_tissue), nrow(aggregated), uniqueN(by_tissue$tissue),
            uniqueN(by_tissue$gene), nrow(background))
)
fwrite(audit, file.path(ref$processed_dir, "reference_audit.tsv"), sep = "\t")

manifest <- data.table(
  resource = c("GTEx v8 CAVIAR", "1000 Genomes Phase3 EUR", "gene mapping", "TSS"),
  source_path = c(ref$caviar_file, ref$background, ref$ensembl_map, ref$tss),
  build = c("hg38", "hg38 lifted from hg19", "Ensembl to HGNC symbol", "hg38"),
  filter = c("PIP > 0.2; |variant-TSS| <= 500kb",
             "EUR MAF > 0.01", "strip Ensembl version then exact map", "one TSS per symbol")
)
fwrite(manifest, file.path(ref$processed_dir, "reference_manifest.tsv"), sep = "\t")
message("Prepared ", uniqueN(by_tissue$tissue), " tissues; ", nrow(by_tissue),
        " tissue-specific and ", nrow(aggregated), " aggregated variant-gene rows")

