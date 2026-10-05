#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
})
source("scripts/eQTL_GTEx/gtex_helpers.R")

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args)) args[[1L]] else "config/gtex_eqtl_validation.yml"
only_dataset <- if (length(args) >= 2L && nzchar(args[[2L]])) args[[2L]] else NULL
config <- read_yaml(config_path)
ref <- config$reference
out_root <- "results/validation/eQTL/GTEx/eQTL1"
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

aggregate_eqtl <- as.data.table(readRDS(file.path(
  ref$processed_dir, "GTEx_v8_CAVIAR_PIPgt0.2_49tissues_aggregated.rds")))
tissue_eqtl <- as.data.table(readRDS(file.path(
  ref$processed_dir, "GTEx_v8_CAVIAR_PIPgt0.2_by_tissue.rds")))
background <- background_for_eqtl(ref$background)
tss <- as.data.table(readRDS(ref$tss))
method_specs <- gtex_method_specs()
n_boot <- as.integer(ref$bootstrap_replicates)
seed <- as.integer(ref$bootstrap_seed)

aggregate_total <- aggregate_eqtl[, .(causal_total = uniqueN(variant_id)), by = gene]
tissue_total <- tissue_eqtl[, .(causal_total = uniqueN(variant_id)), by = .(tissue, gene)]

all_summary <- list()
all_tissue_summary <- list()
audit_rows <- list()

run_scope <- function(dataset, celltype = NULL) {
  scope <- if (is.null(celltype)) "dataset_pool" else "panel"
  scope_label <- if (is.null(celltype)) dataset else paste(dataset, celltype, sep = " / ")
  message("eQTL1 ", scope, ": ", scope_label)
  universe <- read_scope_universe(dataset, config, celltype)
  scope_dir <- if (is.null(celltype)) file.path(out_root, "dataset_pool", dataset) else
    file.path(out_root, "panel", dataset, celltype)
  cache <- scope_cache(
    universe, aggregate_eqtl, tissue_eqtl, background,
    file.path(scope_dir, "_cache"), tss, as.integer(ref$cis_window_bp)
  )
  genes <- unique(universe$gene)

  for (i in seq_len(nrow(method_specs))) {
    spec <- method_specs[i]
    method_dir <- file.path(scope_dir, spec$method)
    complete_files <- file.path(method_dir, c(
      "per_gene_counts.tsv.gz", "per_gene_scent_enrichment.tsv.gz",
      "gene_bootstrap_pooled_or_and_scent_mean.tsv.gz", "summary.tsv",
      "tissue_per_gene_counts.tsv.gz", "tissue_summary.tsv.gz"
    ))
    if (all(file.exists(complete_files)) && all(file.info(complete_files)$size > 0)) {
      message("  ", spec$method, ": existing complete result; skip")
      all_summary[[paste(scope, dataset, celltype, spec$method, sep = "::")]] <<-
        fread(file.path(method_dir, "summary.tsv"))
      all_tissue_summary[[paste(scope, dataset, celltype, spec$method, sep = "::")]] <<-
        fread(file.path(method_dir, "tissue_summary.tsv.gz"))
      next
    }
    selected <- read_scope_selected(dataset, spec, config, celltype)
    if (!nrow(selected)) {
      message("  SKIP ", spec$method, ": no significant links")
      next
    }
    message("  ", spec$method, ": ", nrow(selected), " unique selected links")
    bg_hits <- background_hit_counts(cache$background_hits, selected)
    agg_hits <- selected_hit_counts(cache$aggregate_hits, selected, FALSE)
    counts <- assemble_counts(aggregate_total, agg_hits, bg_hits,
                              cache$background_cis, genes, FALSE)
    result <- summarize_scent_counts(counts, n_boot, seed + i)
    dir.create(method_dir, recursive = TRUE, showWarnings = FALSE)
    fwrite(counts, file.path(method_dir, "per_gene_counts.tsv.gz"), sep = "\t")
    fwrite(result$per_gene,
           file.path(method_dir, "per_gene_scent_enrichment.tsv.gz"), sep = "\t")
    fwrite(result$bootstrap,
           file.path(method_dir, "gene_bootstrap_pooled_or_and_scent_mean.tsv.gz"), sep = "\t")

    summary <- copy(result$summary)
    summary[, `:=`(scope = scope, dataset = dataset,
                   celltype = if (is.null(celltype)) "POOL" else celltype,
                   method = spec$method, n_selected_pairs = nrow(selected),
                   reference = "GTEx_v8_CAVIAR_49_tissues_aggregated_PIPgt0.2",
                   bootstrap_replicates = n_boot)]
    setcolorder(summary, c("scope", "dataset", "celltype", "method"))
    fwrite(summary, file.path(method_dir, "summary.tsv"), sep = "\t")
    all_summary[[paste(scope, dataset, celltype, spec$method, sep = "::")]] <<- summary

    tissue_hits <- selected_hit_counts(cache$tissue_hits, selected, TRUE)
    tissue_counts <- assemble_counts(tissue_total, tissue_hits, bg_hits,
                                     cache$background_cis, genes, TRUE)
    tissue_counts[, `:=`(
      causal_fraction = causal_annotated / causal_total,
      common_fraction = common_annotated / common_total,
      scent_enrichment = fifelse(common_annotated > 0,
                                 (causal_annotated / causal_total) /
                                   (common_annotated / common_total), NA_real_)
    )]
    tissue_summary <- tissue_counts[, .(
      n_genes_all = .N,
      n_genes_scent = sum(is.finite(scent_enrichment)),
      scent_mean = mean(scent_enrichment[is.finite(scent_enrichment)]),
      scent_median = median(scent_enrichment[is.finite(scent_enrichment)]),
      causal_total = sum(causal_total), causal_annotated = sum(causal_annotated),
      common_total = sum(common_total), common_annotated = sum(common_annotated)
    ), by = tissue]
    tissue_summary[, pooled_ratio :=
      (causal_annotated / pmax(common_annotated, 1L)) /
      (causal_total / common_total)]
    tissue_summary[, `:=`(scope = scope, dataset = dataset,
                          celltype = if (is.null(celltype)) "POOL" else celltype,
                          method = spec$method)]
    setcolorder(tissue_summary, c("scope", "dataset", "celltype", "method", "tissue"))
    fwrite(tissue_counts, file.path(method_dir, "tissue_per_gene_counts.tsv.gz"), sep = "\t")
    fwrite(tissue_summary, file.path(method_dir, "tissue_summary.tsv.gz"), sep = "\t")
    all_tissue_summary[[paste(scope, dataset, celltype, spec$method, sep = "::")]] <<-
      tissue_summary
  }
  scope_audit <- data.table(
    scope = scope, dataset = dataset,
    celltype = if (is.null(celltype)) "POOL" else celltype,
    n_candidate_pairs = nrow(universe), n_candidate_genes = uniqueN(universe$gene),
    n_candidate_peaks = uniqueN(universe$peak),
    n_background_peak_overlaps = nrow(cache$background_hits),
    n_aggregate_eqtl_peak_overlaps = nrow(cache$aggregate_hits),
    n_tissue_eqtl_peak_overlaps = nrow(cache$tissue_hits)
  )
  fwrite(scope_audit, file.path(scope_dir, "scope_audit.tsv"), sep = "\t")
  audit_rows[[paste(scope, dataset, celltype, sep = "::")]] <<- scope_audit
}

panels <- gtex_panels(config)
if (!is.null(only_dataset)) panels <- panels[dataset == only_dataset]
for (i in seq_len(nrow(panels))) run_scope(panels$dataset[[i]], panels$celltype[[i]])
datasets_to_run <- names(config$benchmarks)
if (!is.null(only_dataset)) datasets_to_run <- intersect(datasets_to_run, only_dataset)
for (dataset in datasets_to_run) run_scope(dataset, NULL)

summary <- rbindlist(all_summary, fill = TRUE)
summary[, pooled_or_padj := p.adjust(pooled_or_pval, method = "BH")]
suffix <- if (is.null(only_dataset)) "all" else paste0("partial_", only_dataset)
fwrite(summary, file.path(out_root, paste0(suffix, "_results.tsv")), sep = "\t")
fwrite(rbindlist(all_tissue_summary, fill = TRUE),
       file.path(out_root, paste0(suffix, "_tissue_results.tsv.gz")), sep = "\t")
fwrite(rbindlist(audit_rows, fill = TRUE),
       file.path(out_root, paste0(suffix, "_scope_audit.tsv")), sep = "\t")
message("GTEx eQTL1 complete: ", out_root)
