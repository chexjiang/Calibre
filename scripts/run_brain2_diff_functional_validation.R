#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(readxl)
  library(DBI)
  library(RSQLite)
})

args <- commandArgs(trailingOnly = TRUE)
project_dir <- normalizePath(if (length(args) >= 1L) args[[1L]] else ".")
output_dir <- if (length(args) >= 2L) args[[2L]] else
  file.path(project_dir, "results", "validation", "functional", "Brain2_diff")

scmultimap_file <- file.path(
  project_dir, "results", "real", "Brain2_diff", "Microglia",
  "scmultimap", "all_pairs.tsv.gz"
)
epilink_file <- file.path(
  project_dir, "results", "real", "Brain2_diff", "Microglia",
  "epilink_diff_poisson_no_donor_effect_no_donor_robust", "all_pairs.tsv.gz"
)
anderson_file <- file.path(
  project_dir, "reference", "functional", "Brain2_diff",
  "anderson_2023_table_s3.xlsx"
)
mathys_file <- file.path(
  project_dir, "reference", "functional", "Brain2_diff",
  "mathys_2024_microglia_differential_expression.tsv.gz"
)

input_files <- c(scmultimap_file, epilink_file, anderson_file, mathys_file)
missing_files <- input_files[!file.exists(input_files)]
if (length(missing_files)) {
  stop("Missing input file(s):\n", paste(missing_files, collapse = "\n"))
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
figure_dir <- file.path(output_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

clean_gene <- function(x) toupper(trimws(as.character(x)))

scmultimap <- fread(scmultimap_file)
epilink <- fread(epilink_file)
scmultimap[, pair_id := paste(gene, peak, sep = "|")]
epilink[, pair_id := paste(gene, peak, sep = "|")]

common_pair_ids <- intersect(scmultimap$pair_id, epilink$pair_id)
common_genes <- sort(unique(clean_gene(
  scmultimap[pair_id %chin% common_pair_ids, gene]
)))

selected_pairs <- list(
  scMultiMap_padj_0.10 = scmultimap[
    pair_id %chin% common_pair_ids & !is.na(padj) & padj <= 0.10
  ],
  scMultiMap_raw_p_0.05 = scmultimap[
    pair_id %chin% common_pair_ids & !is.na(pval) & pval < 0.05
  ],
  EpiLink_significant = epilink[
    pair_id %chin% common_pair_ids & !is.na(significant) & significant
  ]
)
selected_genes <- lapply(selected_pairs, function(x) {
  sort(unique(clean_gene(x$gene)))
})

for (method in names(selected_pairs)) {
  fwrite(
    selected_pairs[[method]],
    file.path(output_dir, paste0(method, "_selected_pairs.tsv.gz")),
    sep = "\t"
  )
  fwrite(
    data.table(gene = selected_genes[[method]]),
    file.path(output_dir, paste0(method, "_selected_genes.tsv")),
    sep = "\t"
  )
}

pair_summary <- rbindlist(lapply(names(selected_pairs), function(method) {
  data.table(
    method = method,
    tested_pairs = length(common_pair_ids),
    selected_pairs = nrow(selected_pairs[[method]]),
    selected_pair_fraction = nrow(selected_pairs[[method]]) / length(common_pair_ids),
    tested_genes = length(common_genes),
    selected_genes = length(selected_genes[[method]])
  )
}))
fwrite(pair_summary, file.path(output_dir, "pair_and_gene_summary.tsv"), sep = "\t")

pair_overlap <- rbindlist(lapply(combn(names(selected_pairs), 2L, simplify = FALSE), function(pair) {
  a <- selected_pairs[[pair[[1L]]]]$pair_id
  b <- selected_pairs[[pair[[2L]]]]$pair_id
  shared <- length(intersect(a, b))
  data.table(
    method_a = pair[[1L]],
    method_b = pair[[2L]],
    pairs_a = length(a),
    pairs_b = length(b),
    overlap_pairs = shared,
    jaccard = shared / length(union(a, b))
  )
}))
fwrite(pair_overlap, file.path(output_dir, "pair_overlap.tsv"), sep = "\t")

summary_long <- melt(
  pair_summary,
  id.vars = "method",
  measure.vars = c("selected_pairs", "selected_genes"),
  variable.name = "quantity",
  value.name = "count"
)
summary_long[, method_label := fcase(
  method == "scMultiMap_padj_0.10", "scMultiMap (padj <= 0.10)",
  method == "scMultiMap_raw_p_0.05", "scMultiMap (raw p < 0.05)",
  default = "EpiLink-Diff Poisson"
)]
p_counts <- ggplot(summary_long, aes(method_label, count, fill = quantity)) +
  geom_col(position = "dodge") +
  labs(x = NULL, y = "Count", fill = NULL) +
  facet_wrap(~quantity, scales = "free_y") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "top")
ggsave(file.path(figure_dir, "differential_pair_gene_counts.png"), p_counts,
       width = 8.5, height = 4.2, dpi = 300)
ggsave(file.path(figure_dir, "differential_pair_gene_counts.pdf"), p_counts,
       width = 8.5, height = 4.2)

anderson <- as.data.table(read_excel(anderson_file, sheet = "TableS3_ADCtrl_DEGs"))
anderson_degs <- sort(unique(clean_gene(anderson[
  celltype == "Microglia" & p_val_adj < 0.01 & abs(avg_log2FC) > 0.25,
  gene
])))

mathys <- fread(mathys_file)
mathys_degs <- sort(unique(clean_gene(mathys[
  path == "nrad" & region == "allregions" & col_nm != 0,
  gene
])))

deg_sets <- list(
  Anderson_2023_Microglia = anderson_degs,
  Mathys_2024_Microglia_NIA_Reagan_all_regions = mathys_degs
)

fwrite(
  rbindlist(lapply(names(deg_sets), function(reference) {
    data.table(reference = reference, gene = deg_sets[[reference]])
  })),
  file.path(output_dir, "reference_deg_sets.tsv.gz"),
  sep = "\t"
)

membership <- data.table(gene = common_genes)
for (method in names(selected_genes)) {
  membership[, (method) := gene %chin% selected_genes[[method]]]
}
for (reference in names(deg_sets)) {
  membership[, (reference) := gene %chin% deg_sets[[reference]]]
}
fwrite(
  membership,
  file.path(output_dir, "tested_gene_membership.tsv.gz"),
  sep = "\t"
)

enrichment_rows <- list()
overlap_rows <- list()
k <- 1L
for (method in names(selected_genes)) {
  method_genes <- intersect(selected_genes[[method]], common_genes)
  for (reference in names(deg_sets)) {
    reference_genes <- intersect(deg_sets[[reference]], common_genes)
    overlap <- intersect(method_genes, reference_genes)
    tab <- matrix(c(
      length(overlap),
      length(setdiff(method_genes, reference_genes)),
      length(setdiff(reference_genes, method_genes)),
      length(setdiff(common_genes, union(method_genes, reference_genes)))
    ), nrow = 2L, byrow = TRUE)
    fit <- fisher.test(tab, alternative = "greater")
    enrichment_rows[[k]] <- data.table(
      method = method,
      reference = reference,
      universe_genes = length(common_genes),
      selected_genes = length(method_genes),
      reference_degs_in_universe = length(reference_genes),
      overlap_genes = length(overlap),
      odds_ratio = unname(fit$estimate),
      conf_low = fit$conf.int[[1L]],
      conf_high = fit$conf.int[[2L]],
      fisher_p = fit$p.value
    )
    overlap_rows[[k]] <- data.table(
      method = method,
      reference = reference,
      gene = overlap
    )
    k <- k + 1L
  }
}
deg_enrichment <- rbindlist(enrichment_rows)
deg_enrichment[, fisher_fdr := p.adjust(fisher_p, method = "BH")]
fwrite(deg_enrichment, file.path(output_dir, "deg_enrichment.tsv"), sep = "\t")
fwrite(rbindlist(overlap_rows), file.path(output_dir, "deg_overlap_genes.tsv"), sep = "\t")

plot_enrichment <- copy(deg_enrichment)
plot_enrichment[, method_label := fcase(
  method == "scMultiMap_padj_0.10", "scMultiMap (padj <= 0.10)",
  method == "scMultiMap_raw_p_0.05", "scMultiMap (raw p < 0.05)",
  default = "EpiLink-Diff Poisson"
)]
plot_enrichment[, method_label := factor(
  method_label,
  levels = c(
    "EpiLink-Diff Poisson",
    "scMultiMap (raw p < 0.05)",
    "scMultiMap (padj <= 0.10)"
  )
)]
plot_enrichment[, reference := sub("_Microglia.*$", "", reference)]
p_deg <- ggplot(plot_enrichment, aes(odds_ratio, method_label, colour = reference)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey55") +
  geom_errorbar(aes(xmin = conf_low, xmax = conf_high), width = 0.16,
                orientation = "y", position = position_dodge(width = 0.45)) +
  geom_point(size = 2.6, position = position_dodge(width = 0.45)) +
  scale_x_log10() +
  labs(x = "Odds ratio (log scale)", y = NULL, colour = "DEG reference") +
  theme_bw(base_size = 11) +
  theme(legend.position = "top")
ggsave(file.path(figure_dir, "deg_enrichment.pdf"), p_deg, width = 7.3, height = 3.2)
ggsave(file.path(figure_dir, "deg_enrichment.png"), p_deg, width = 7.3, height = 3.2,
       dpi = 300)

message("Loading local human GO:BP annotations")
go_con <- dbConnect(
  SQLite(),
  system.file("extdata", "GO.sqlite", package = "GO.db")
)
go_term_names <- as.data.table(dbGetQuery(
  go_con,
  "SELECT go_id AS term_id, term AS term_name FROM go_term WHERE ontology = 'BP'"
))
dbDisconnect(go_con)

go_enrichment <- function(query, background = NULL) {
  if (is.null(background)) stop("A tested-gene background is required")
  background <- sort(unique(clean_gene(background)))
  query <- intersect(sort(unique(clean_gene(query))), background)
  if (length(query) < 2L) return(data.table())

  con <- dbConnect(
    SQLite(),
    system.file("extdata", "org.Hs.eg.sqlite", package = "org.Hs.eg.db")
  )
  on.exit(dbDisconnect(con), add = TRUE)
  dbWriteTable(
    con,
    "analysis_background",
    data.frame(gene = background, is_query = as.integer(background %in% query)),
    temporary = TRUE,
    overwrite = TRUE
  )
  dbExecute(con, "CREATE INDEX background_gene_idx ON analysis_background(gene)")

  sizes <- dbGetQuery(con, paste(
    "SELECT",
    "COUNT(DISTINCT b.gene) AS universe_size,",
    "COUNT(DISTINCT CASE WHEN b.is_query = 1 THEN b.gene END) AS query_size",
    "FROM analysis_background b",
    "JOIN gene_info gi ON UPPER(gi.symbol) = b.gene",
    "JOIN go_bp_all gb ON gb._id = gi._id"
  ))
  result <- as.data.table(dbGetQuery(con, paste(
    "SELECT gb.go_id AS term_id,",
    "COUNT(DISTINCT b.gene) AS term_size,",
    "COUNT(DISTINCT CASE WHEN b.is_query = 1 THEN b.gene END) AS intersection_size,",
    "GROUP_CONCAT(DISTINCT CASE WHEN b.is_query = 1 THEN b.gene END) AS intersection_genes",
    "FROM analysis_background b",
    "JOIN gene_info gi ON UPPER(gi.symbol) = b.gene",
    "JOIN go_bp_all gb ON gb._id = gi._id",
    "GROUP BY gb.go_id",
    "HAVING intersection_size > 0"
  )))
  if (!nrow(result)) return(result)
  result[, query_size := sizes$query_size[[1L]]]
  result[, precision := intersection_size / query_size]
  result[, recall := intersection_size / term_size]
  result[, p_value := phyper(
    intersection_size - 1L,
    term_size,
    sizes$universe_size[[1L]] - term_size,
    query_size,
    lower.tail = FALSE
  )]
  result[, fdr := p.adjust(p_value, method = "BH")]
  result[, significant := fdr < 0.05]
  result <- merge(result, go_term_names, by = "term_id", all.x = TRUE, sort = FALSE)
  result[order(p_value)]
}

go_query_rows <- rbindlist(lapply(names(selected_genes), function(method) {
  rbind(
    data.table(
      method = method,
      analysis = "direct_selected_genes",
      gene = selected_genes[[method]]
    ),
    data.table(
      method = method,
      analysis = "article_like_anderson_intersection",
      gene = intersect(selected_genes[[method]], anderson_degs)
    )
  )
}))
fwrite(
  go_query_rows,
  file.path(output_dir, "go_input_genes_direct_and_article_like.tsv"),
  sep = "\t"
)

go_results <- list()
for (analysis_name in unique(go_query_rows$analysis)) {
  for (method_name in names(selected_genes)) {
    query <- go_query_rows[
      method == method_name & analysis == analysis_name, gene
    ]
    message("Local GO:BP: ", method_name, " [", analysis_name, "]")
    ans <- go_enrichment(query, common_genes)
    if (nrow(ans)) {
      ans[, `:=`(method = method_name, analysis = analysis_name,
                 background = "tested_candidate_genes")]
      go_results[[paste(method_name, analysis_name, sep = "|")]] <- ans
    }
  }
}
go_results <- rbindlist(go_results, fill = TRUE)
fwrite(go_results, file.path(output_dir, "go_bp_enrichment.tsv.gz"), sep = "\t")
fwrite(
  go_results[analysis == "direct_selected_genes"],
  file.path(output_dir, "direct_go_enrichment.tsv.gz"),
  sep = "\t"
)
fwrite(
  go_results[analysis == "article_like_anderson_intersection"],
  file.path(output_dir, "article_like_go_enrichment.tsv.gz"),
  sep = "\t"
)

plot_go_results <- function(go_data, stem) {
  go_plot_data <- go_data[significant == TRUE]
  if (!nrow(go_plot_data)) return(invisible(NULL))
  go_plot_data[, neg_log10_fdr := -log10(fdr)]
  go_plot_data <- go_plot_data[order(method, fdr), head(.SD, 10L), by = method]
  go_plot_data[, method_label := fcase(
    method == "scMultiMap_padj_0.10", "scMultiMap (padj)",
    method == "scMultiMap_raw_p_0.05", "scMultiMap (raw p)",
    default = "EpiLink-Diff Poisson"
  )]
  go_plot_data[, plot_label := paste0(term_name, " (", term_id, ")")]
  go_plot_data[, plot_label := factor(
    plot_label, levels = rev(unique(plot_label[order(neg_log10_fdr)]))
  )]
  p_go <- ggplot(go_plot_data, aes(neg_log10_fdr, plot_label, fill = method_label)) +
    geom_col(show.legend = FALSE) +
    facet_wrap(~method_label, scales = "free_y") +
    labs(x = expression(-log[10]("BH FDR")), y = NULL) +
    theme_bw(base_size = 10)
  ggsave(file.path(figure_dir, paste0(stem, ".pdf")), p_go, width = 10, height = 6)
  ggsave(file.path(figure_dir, paste0(stem, ".png")), p_go, width = 10, height = 6,
         dpi = 300)
}
plot_go_results(
  go_results[analysis == "direct_selected_genes"],
  "direct_go_enrichment"
)
plot_go_results(
  go_results[analysis == "article_like_anderson_intersection"],
  "article_like_go_enrichment"
)

manifest <- data.table(
  item = c(
    "scMultiMap padj selection", "scMultiMap raw-p selection", "EpiLink selection", "gene universe",
    "Anderson DEG definition", "Mathys DEG definition", "DEG enrichment test",
    "direct GO query", "article-like GO query", "GO correction", "GO background"
  ),
  value = c(
    "padj <= 0.10",
    "raw p < 0.05",
    "significant == TRUE (EpiLink calibrated_q <= configured FDR)",
    "genes in peak-gene pairs tested by both methods",
    "Microglia; p_val_adj < 0.01; abs(avg_log2FC) > 0.25",
    "Mic_Immune_Mic; NIA-Reagan AD status (nrad); allregions; col_nm != 0",
    "one-sided Fisher exact test; BH across six method-reference tests",
    "all selected genes",
    "selected genes intersected with Anderson 2023 Microglia DEGs",
    paste0(
      "local org.Hs.eg.db/GO.db GO:BP; one-sided hypergeometric test; BH FDR; ",
      "org.Hs.eg.db ", as.character(packageVersion("org.Hs.eg.db")), ", GO.db ",
      as.character(packageVersion("GO.db"))
    ),
    "common genes in peak-gene pairs tested by all methods"
  )
)
fwrite(manifest, file.path(output_dir, "analysis_manifest.tsv"), sep = "\t")

readme <- c(
  "# Brain2 differential peak-gene functional validation",
  "",
  "The disease/control differential panels in the scMultiMap paper are Fig. 4b-c;",
  "the literal Fig. 3b-c are simulation/PLAC-seq consistency panels and require",
  "external simulation truth or PLAC-seq reference data. This run therefore produces",
  "the disease/control Fig. 4b-c-like pair/gene comparison and functional analyses.",
  "",
  "Three primary sets are compared: scMultiMap `padj <= 0.10`, scMultiMap raw `p < 0.05`,",
  "and EpiLink-Diff Poisson `significant == TRUE`. The scMultiMap raw set is the",
  "nominal threshold used in the paper; it is not FDR-controlled.",
  "",
  "DEG enrichment is evaluated at the gene level with a one-sided Fisher exact test.",
  "The universe is the set of genes in candidate pairs tested by both methods.",
  "Direct GO:BP uses all selected genes. Article-like GO:BP uses selected genes intersected",
  "with Anderson et al. Microglia DEGs, matching the paper's differential functional analysis.",
  "The GO background is the common tested candidate-gene universe.",
  "GO enrichment uses local Bioconductor annotations; no result-derived genes are uploaded.",
  "",
  paste0("Generated: ", format(Sys.time(), tz = "UTC", usetz = TRUE))
)
writeLines(readme, file.path(output_dir, "README.md"))

message("Results written to: ", normalizePath(output_dir))
