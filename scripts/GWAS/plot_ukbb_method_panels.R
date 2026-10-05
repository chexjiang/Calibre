#!/usr/bin/env Rscript

# UKBB-only method comparison. Each facet is a dataset (pooled mode) or a
# dataset-cell-type (celltype mode); x-axis is method.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly=TRUE)
root <- if (length(args)) normalizePath(args[[1]], mustWork=TRUE) else normalizePath(".", mustWork=TRUE)
mode <- if (length(args) >= 2L) args[[2]] else "pooled"
if (!mode %in% c("pooled", "celltype")) stop("mode must be pooled or celltype")
analysis_name <- if (length(args) >= 3L) args[[3]] else
  if (mode == "pooled") "method_comparison" else "celltype_comparison"

if (mode == "pooled") {
  base <- file.path(root, "results", "validation", "GWAS", "SCENT", "sakaue_style", analysis_name)
  infile <- file.path(base, "results", "scent_sakaue_summary.tsv")
  outfile <- file.path(base, "figures", "scent_sakaue_ukbb_method_by_dataset.pdf")
  tablefile <- file.path(base, "figures", "scent_sakaue_ukbb_method_by_dataset.tsv")
} else {
  base <- file.path(root, "results", "validation", "GWAS", "SCENT", "sakaue_style", analysis_name)
  infile <- file.path(base, "results", "scent_sakaue_summary.tsv")
  outfile <- file.path(base, "figures", "scent_sakaue_ukbb_method_by_celltype.pdf")
  tablefile <- file.path(base, "figures", "scent_sakaue_ukbb_method_by_celltype.tsv")
}
dir.create(dirname(outfile), recursive=TRUE, showWarnings=FALSE)

x <- fread(infile)[source == "UKBB"]
x <- x[is.finite(mean_enrichment) & is.finite(ci_low) & is.finite(ci_high)]
method_levels <- c("EpiLink-NB", "EpiLink-Poisson", "scMultiMap",
                   "Signac (corrected)", "Signac (uncorrected)")
x[, method_label := factor(fifelse(grepl("EpiLink-NB", method_label), "EpiLink-NB",
                                    fifelse(grepl("EpiLink-Poisson", method_label), "EpiLink-Poisson", method_label)),
                            levels=method_levels)]
x[, dataset_label := fifelse(dataset == "Dogma_CTRL", "Dogma CTRL",
                             fifelse(dataset == "Dogma_STIM", "Dogma STIM", dataset))]
if (mode == "pooled") {
  x[, panel_label := factor(dataset_label,
                            levels=c("pbmc", "bmmc", "Dogma CTRL", "Dogma STIM",
                                     "NEAT", "Arthritis", "Brain"))]
  ncol <- 3L; width <- 13; height <- 9
  title <- "UKBB GWAS enrichment by dataset"
  subtitle <- "Pooled over primary cell types; point = mean across 94 traits; whisker = 95% trait-bootstrap CI"
} else {
  x[, panel_label := factor(paste(dataset_label, celltypes, sep="\n"))]
  ncol <- 4L; width <- 16; height <- 15
  title <- "UKBB GWAS enrichment by dataset and cell type"
  subtitle <- "Point = mean across 94 traits; whisker = 95% trait-bootstrap CI"
}

p <- ggplot(x, aes(x=method_label, y=mean_enrichment, color=method_label, shape=method_label)) +
  geom_hline(yintercept=1, linewidth=0.35, linetype="dashed", color="grey45") +
  geom_errorbar(aes(ymin=pmax(ci_low,0), ymax=ci_high), width=0.16, linewidth=0.55) +
  geom_point(size=2.5) +
  facet_wrap(~panel_label, ncol=ncol, scales="free_y", drop=FALSE) +
  scale_color_manual(values=c("EpiLink-NB"="#1b6ca8", "EpiLink-Poisson"="#d97706",
                              "scMultiMap"="#238b45", "Signac (corrected)"="#7b3294",
                              "Signac (uncorrected)"="#c51b7d")) +
  scale_shape_manual(values=c("EpiLink-NB"=16, "EpiLink-Poisson"=17,
                              "scMultiMap"=15, "Signac (corrected)"=18,
                              "Signac (uncorrected)"=8)) +
  scale_y_continuous(expand=expansion(mult=c(0,0.10))) +
  labs(title=title, subtitle=subtitle, x="Method", y="Mean causal-variant enrichment",
       color=NULL, shape=NULL) +
  theme_classic(base_size=11) +
  theme(axis.text.x=element_text(angle=38, hjust=1, vjust=1),
        strip.background=element_blank(), strip.text=element_text(face="bold", size=9),
        legend.position="bottom", plot.title=element_text(face="bold"),
        panel.spacing=unit(0.9,"lines"))

ggsave(outfile, p, width=width, height=height, units="in", device=cairo_pdf)
setorder(x, panel_label, method_label)
fwrite(x, tablefile, sep="\t")
message("Wrote ", outfile)
