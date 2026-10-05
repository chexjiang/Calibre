#!/usr/bin/env Rscript

# Fig. 4-style locus panels for selected EpiLink-only GWAS mechanisms.
# ATAC is shown as pseudo-bulk accessibility over called peaks because the
# processed matrices do not retain base-pair fragment coverage.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(Matrix)
  library(Seurat)
  library(Signac)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[[1L]] else ".", mustWork = TRUE)
out <- file.path(root, "results", "validation", "functional", "GWAS_loci")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

loci <- list(
  IL2RA = list(
    id = "immune_IL2RA_rs61839660", dataset = "Arthritis", celltype = "T",
    chr = "chr10", x0 = 5750000L, x1 = 6350000L, pos = 6052734L,
    variant = "rs61839660", gene = "IL2RA", peak = "chr10-6052701-6052901",
    epilink_dir = "epilink_poisson", model = "EpiLink-Poisson",
    trait = "Inflammatory bowel disease", fine_map = "IBD PIP = 0.999",
    gwas_type = "tabix", gwas_file = file.path(root, "results", "validation", "GWAS",
      "sLDSC", "sumstats", "raw", "finngen_R12_K11_IBD_STRICT.gz")),
  CLU = list(
    id = "brain_CLU_rs1532278", dataset = "Brain2", celltype = "Astrocytes",
    chr = "chr8", x0 = 27200000L, x1 = 28000000L, pos = 27608798L,
    variant = "rs1532278", gene = "CLU", peak = "chr8-27608401-27609282",
    epilink_dir = "epilink_poisson_condition_effect",
    model = "EpiLink-Poisson + condition", trait = "Alzheimer disease",
    fine_map = "FinnGen dementia PIP = 0.296",
    gwas_type = "tabix", gwas_file = file.path(root, "results", "validation", "GWAS",
      "sLDSC", "sumstats", "raw", "finngen_R12_G6_ALZHEIMER.gz")),
  FSHB = list(
    id = "pituitary_FSHB_rs11031006", dataset = "Pituitary", celltype = "Gonadotropes",
    chr = "chr11", x0 = 29700000L, x1 = 30700000L, pos = 30204809L,
    variant = "rs11031006", gene = "FSHB", peak = "chr11-30203591-30205849",
    epilink_dir = "epilink_poisson", model = "EpiLink-Poisson",
    trait = "Endometriosis", fine_map = "FinnGen PIP = 0.520",
    gwas_type = "finemap", gwas_file = file.path(root, "results", "validation", "GWAS",
      "reference", "FinnGen", "raw", "summary",
      "finngen_R12_N14_ENDOMETRIOSIS.SUSIE.snp.filter.tsv"))
  ,
  GH1 = list(
    id = "pituitary_GH1_rs11568828", dataset = "Pituitary", celltype = "Thyrotropes",
    chr = "chr17", x0 = 63500000L, x1 = 64300000L, pos = 63918913L,
    variant = "rs11568828", gene = "GH1", peak = "chr17-63918218-63919459",
    epilink_dir = "epilink_poisson", model = "EpiLink-Poisson",
    trait = "Height", fine_map = "UKBB height causal variant",
    gwas_type = "finemap", gwas_file = file.path(root, "results", "validation", "GWAS",
      "reference", "FinnGen", "raw", "summary",
      "finngen_R12_HEIGHT_IRN.SUSIE.snp.filter.tsv"))
)

parse_peaks <- function(x) {
  z <- tstrsplit(as.character(x), "-", fixed = TRUE)
  data.table(chr = z[[1L]], start = as.integer(z[[2L]]), end = as.integer(z[[3L]]))
}

result_root <- function(cfg) {
  p <- file.path(root, "results", "real", cfg$dataset, cfg$celltype)
  if (cfg$dataset == "Brain2") p <- file.path(p, "ALL")
  p
}

read_pairs <- function(path) {
  if (!file.exists(path)) return(data.table(gene = character(), peak = character()))
  x <- tryCatch(fread(path, showProgress = FALSE), error = function(e) data.table())
  if (!all(c("gene", "peak") %chin% names(x))) return(data.table(gene = character(), peak = character()))
  x
}

read_gwas <- function(cfg) {
  if (cfg$gwas_type == "tabix") {
    region <- sprintf("%s:%d-%d", sub("^chr", "", cfg$chr), cfg$x0, cfg$x1)
    txt <- system2("tabix", c(cfg$gwas_file, region), stdout = TRUE)
    if (!length(txt)) return(data.table())
    x <- fread(text = paste(txt, collapse = "\n"), header = FALSE, sep = "\t")
    setnames(x, c("chrom", "pos", "ref", "alt", "rsid", "nearest_gene", "p",
                  "mlogp", "beta", "se", "af", "af_cases", "af_controls"))
    x[, pos := as.integer(pos)]
    x[, highlighted := grepl(paste0("(^|,)", cfg$variant, "($|,)"), rsid)]
  } else {
    x <- fread(cfg$gwas_file, showProgress = FALSE)
    x <- x[chromosome == cfg$chr & position >= cfg$x0 & position <= cfg$x1]
    x[, `:=`(pos = as.integer(position), mlogp = -log10(p), rsid = v,
             highlighted = position == cfg$pos)]
  }
  x[is.finite(mlogp)]
}

extract_object_tracks <- function(cfg) {
  obj_file <- file.path(result_root(cfg), "signac", "signac_object.rds")
  if (!file.exists(obj_file)) stop("Missing ", obj_file)
  obj <- readRDS(obj_file)
  counts <- GetAssayData(obj, assay = "ATAC", layer = "counts")
  pp <- parse_peaks(rownames(counts))
  keep <- which(pp$chr == cfg$chr & pp$end >= cfg$x0 & pp$start <= cfg$x1)
  atac <- pp[keep]
  if (length(keep)) {
    denom <- sum(Matrix::colSums(counts)) / 1e6
    atac[, cpm := Matrix::rowSums(counts[keep, , drop = FALSE]) / denom]
  } else atac[, cpm := numeric()]

  gr <- Signac::Annotation(obj[["ATAC"]])
  take <- as.character(GenomicRanges::seqnames(gr)) == cfg$chr &
    IRanges::end(gr) >= cfg$x0 & IRanges::start(gr) <= cfg$x1
  gr <- gr[take]
  ann <- data.table(chr = as.character(GenomicRanges::seqnames(gr)),
                    start = IRanges::start(gr), end = IRanges::end(gr),
                    strand = as.character(GenomicRanges::strand(gr)), gene = as.character(gr$gene_name),
                    type = as.character(gr$type))
  ann <- ann[!is.na(gene) & nzchar(gene)]
  genes <- ann[, .(start = min(start), end = max(end), strand = strand[1L]), by = gene]
  rm(obj, counts, gr); gc(verbose = FALSE)
  list(atac = atac, genes = genes)
}

method_tracks <- function(cfg) {
  rr <- result_root(cfg)
  paths <- list(
    EpiLink = file.path(rr, cfg$epilink_dir, "all_pairs.tsv.gz"),
    scMultiMap = file.path(rr, "scmultimap", "significant_pairs.tsv.gz"),
    `Signac corrected` = file.path(rr, "signac", "significant_pairs_corrected.tsv.gz"),
    `Signac uncorrected` = file.path(rr, "signac", "significant_pairs_uncorrected.tsv.gz")
  )
  ans <- rbindlist(lapply(names(paths), function(method) {
    x <- read_pairs(paths[[method]])
    if (!nrow(x)) return(NULL)
    if (method == "EpiLink") x <- x[is.finite(padj) & padj < 0.1]
    x <- x[gene == cfg$gene]
    if (!nrow(x)) return(NULL)
    cbind(x[, .(gene, peak)], parse_peaks(x$peak))[, method := method]
  }), fill = TRUE)
  if (!nrow(ans)) ans <- data.table(gene = character(), peak = character(), chr = character(),
                                    start = integer(), end = integer(), method = character())
  ans[chr == cfg$chr & end >= cfg$x0 & start <= cfg$x1]
}

target_stats <- function(cfg) {
  rr <- result_root(cfg)
  paths <- list(
    EpiLink = file.path(rr, cfg$epilink_dir, "all_pairs.tsv.gz"),
    scMultiMap = file.path(rr, "scmultimap", "all_pairs.tsv.gz"),
    `Signac corrected` = file.path(rr, "signac", "all_pairs.tsv.gz"),
    `Signac uncorrected` = file.path(rr, "signac", "all_pairs.tsv.gz")
  )
  ans <- rbindlist(lapply(names(paths), function(method) {
    x <- read_pairs(paths[[method]])
    x <- x[gene == cfg$gene & peak == cfg$peak]
    if (!nrow(x)) return(data.table(method = method, p = NA_real_, padj = NA_real_))
    if (method == "Signac corrected") {
      p_candidates <- c("pvalue_corrected", "pvalue", "pval")
      q_candidates <- c("padj_corrected", "padj")
    } else if (method == "Signac uncorrected") {
      p_candidates <- c("pvalue_uncorrected", "pvalue", "pval")
      q_candidates <- c("padj_uncorrected", "padj")
    } else {
      p_candidates <- c("pval", "pvalue")
      q_candidates <- c("padj")
    }
    pcol <- intersect(p_candidates, names(x))
    qcol <- intersect(q_candidates, names(x))
    data.table(method = method,
               p = if (length(pcol)) as.numeric(x[[pcol[[1L]]]][1L]) else NA_real_,
               padj = if (length(qcol)) as.numeric(x[[qcol[[1L]]]][1L]) else NA_real_)
  }), fill = TRUE)
  ans[, significant := is.finite(padj) & padj < 0.1]
  ans
}

make_arc <- function(a, n = 80L) {
  tt <- seq(0, 1, length.out = n)
  h <- pmin(0.9, pmax(0.18, abs(a$tss - a$mid) / 3e5))
  data.table(link = a$link, x = (1 - tt) * a$mid + tt * a$tss,
             y = -4 * h * tt * (1 - tt), highlighted = a$highlighted)
}

plot_one <- function(cfg) {
  message("Plotting ", cfg$id)
  gwas <- read_gwas(cfg)
  tr <- extract_object_tracks(cfg)
  mt <- method_tracks(cfg)
  stats <- target_stats(cfg)
  stats[, locus_id := cfg$id]
  fwrite(stats, file.path(out, paste0(cfg$id, "_method_stats.tsv")), sep = "\t")

  rr <- result_root(cfg)
  epi <- read_pairs(file.path(rr, cfg$epilink_dir, "all_pairs.tsv.gz"))
  if (nrow(epi)) {
    epi <- epi[is.finite(padj) & padj < 0.1]
    q <- parse_peaks(epi$peak)
    epi <- cbind(epi, q)[chr == cfg$chr & end >= cfg$x0 & start <= cfg$x1]
  }
  genes <- tr$genes
  genes[, tss := ifelse(strand == "+", start, end)]
  arcs <- merge(epi, genes[, .(gene, tss)], by = "gene", all = FALSE)
  arcs[, `:=`(mid = (start + end) / 2, link = .I,
              highlighted = gene == cfg$gene & peak == cfg$peak)]
  arc_df <- if (nrow(arcs)) rbindlist(lapply(seq_len(nrow(arcs)), function(i) make_arc(arcs[i]))) else data.table()

  shade <- annotate("rect", xmin = cfg$pos - 2500, xmax = cfg$pos + 2500,
                    ymin = -Inf, ymax = Inf, fill = "#56B4E9", alpha = 0.10)
  xscale <- scale_x_continuous(limits = c(cfg$x0, cfg$x1),
                               labels = function(x) sprintf("%.2f", x / 1e6))
  base_theme <- theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

  p1 <- ggplot(gwas, aes(pos, mlogp)) + shade +
    geom_point(aes(colour = highlighted), size = 0.75, alpha = 0.75) +
    geom_hline(yintercept = -log10(5e-8), linetype = 3, colour = "#B2182B") +
    geom_vline(xintercept = cfg$pos, colour = "#0072B2", linewidth = 0.45) +
    annotate("text", x = cfg$pos, y = Inf, label = cfg$variant, vjust = 1.2,
             hjust = -0.05, colour = "#0072B2", size = 3) +
    scale_colour_manual(values = c(`FALSE` = "grey45", `TRUE` = "#0072B2"), guide = "none") +
    xscale + labs(title = paste0("GWAS: ", cfg$trait), y = expression(-log[10](italic(P))), x = NULL) + base_theme

  atac <- tr$atac
  p2 <- ggplot(atac) + shade +
    geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = cpm),
              fill = "#4D4D4D", colour = NA) + xscale +
    labs(title = paste0("Pseudo-bulk ATAC peaks: ", cfg$dataset, " ", cfg$celltype),
         y = "CPM", x = NULL) + base_theme

  method_levels <- c("EpiLink", "scMultiMap", "Signac corrected", "Signac uncorrected")
  mt[, method := factor(method, levels = rev(method_levels))]
  p3 <- ggplot(mt) + shade +
    geom_rect(aes(xmin = start, xmax = end,
                  ymin = as.numeric(method) - 0.32, ymax = as.numeric(method) + 0.32,
                  fill = method), colour = NA) +
    xscale + scale_y_continuous(breaks = seq_along(method_levels), labels = rev(method_levels)) +
    scale_fill_manual(values = c(EpiLink = "#009E73", scMultiMap = "#CC79A7",
                                 `Signac corrected` = "#E69F00", `Signac uncorrected` = "#F0E442"),
                      guide = "none", drop = FALSE) +
    labs(title = paste0(cfg$gene, "-linked significant peaks"), y = NULL, x = NULL) + base_theme

  p4 <- ggplot() + shade + geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.25) +
    {if (nrow(arc_df)) geom_path(data = arc_df[highlighted == FALSE],
       aes(x, y, group = link), colour = "#009E73", alpha = 0.30, linewidth = 0.35)} +
    {if (nrow(arc_df)) geom_path(data = arc_df[highlighted == TRUE],
       aes(x, y, group = link), colour = "#0072B2", linewidth = 1.2)} +
    xscale + scale_y_continuous(breaks = NULL) +
    labs(title = paste0(cfg$model, " local links (blue = ", cfg$gene, " target)"), y = NULL, x = NULL) + base_theme

  genes <- genes[end >= cfg$x0 & start <= cfg$x1]
  genes[, y := frank(start, ties.method = "dense") %% 3L]
  p5 <- ggplot(genes) + shade +
    geom_segment(aes(x = start, xend = end, y = y, yend = y), linewidth = 0.5,
                 arrow = arrow(length = unit(0.05, "inches"), type = "closed"), colour = "#2C7FB8") +
    geom_text(aes(x = (start + end) / 2, y = y + 0.28, label = gene), size = 2.5,
              check_overlap = TRUE) + xscale + scale_y_continuous(breaks = NULL) +
    labs(title = "Gene bodies", y = NULL,
         x = paste0(cfg$chr, " position (Mb)")) + base_theme

  fig <- p1 / p2 / p3 / p4 / p5 +
    plot_layout(heights = c(2.2, 1.0, 1.2, 1.15, 1.1)) +
    plot_annotation(
      title = paste0(cfg$variant, " – ", cfg$gene, " (", cfg$fine_map, ")"),
      subtitle = "GRCh38; blue stripe marks the fine-mapped variant. Regional points are not LD-coloured.",
      theme = theme(plot.title = element_text(face = "bold", size = 14),
                    plot.subtitle = element_text(size = 9)))
  ggsave(file.path(out, paste0(cfg$id, ".png")), fig, width = 10.5, height = 10, dpi = 300)
  ggsave(file.path(out, paste0(cfg$id, ".pdf")), fig, width = 10.5, height = 10,
         device = cairo_pdf)
  fwrite(gwas, file.path(out, paste0(cfg$id, "_gwas.tsv.gz")), sep = "\t")
  fwrite(mt, file.path(out, paste0(cfg$id, "_method_tracks.tsv.gz")), sep = "\t")
}

for (cfg in loci) plot_one(cfg)
message("Wrote locus panels to ", out)
