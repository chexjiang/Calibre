#!/usr/bin/env Rscript

suppressPackageStartupMessages({library(data.table); library(ggplot2)})
args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(if (length(args)) args[[1]] else ".", mustWork=TRUE)
source(file.path(root, "scripts", "GWAS", "sLDSC", "source_results_paths.R"))
paths <- sldsc_source_results_paths(root)
dir.create(paths$figures, recursive=TRUE, showWarnings=FALSE)
dir.create(paths$tables, recursive=TRUE, showWarnings=FALSE)

read_view <- function(view_name) {
  f <- file.path(paths[[view_name]], "results", "sldsc_all_results.tsv")
  if (!file.exists(f)) stop("Missing ", f)
  x <- fread(f)
  x[, view := view_name]
  x
}
x <- rbindlist(list(read_view("pooled"), read_view("celltype")), fill=TRUE)

method_order <- c("EpiLink-Poisson", "EpiLink-NB", "scMultiMap",
                  "Signac (corrected)", "Signac (uncorrected)")
method_colors <- c(
  "EpiLink-Poisson"="#E47159", "EpiLink-NB"="#C87878",
  "scMultiMap"="#edd089", "Signac (corrected)"="#54b2a2",
  "Signac (uncorrected)"="#78A878"
)
x[, method_label := factor(method_label, levels=method_order)]
x[, enrichment_low := Enrichment - Enrichment_std_error]
x[, enrichment_high := Enrichment + Enrichment_std_error]
x[, significance := fifelse(Enrichment_p < 0.01, "**", fifelse(Enrichment_p < 0.1, "*", ""))]
x[, matched_trait :=
    (panel == "brain" & trait == "G6_ALZHEIMER") |
    (panel == "immune" & trait %chin% c("M13_RHEUMA", "T1D", "K11_IBD_STRICT", "G6_MS")) |
    (panel == "arthritis_context" & trait == "M13_RHEUMA")]

fwrite(x, file.path(paths$tables, "sldsc_all_results.tsv"), sep="\t")
counts <- x[, .(
  n_results=.N, n_annotations=uniqueN(annotation_id),
  n_significant_p01=sum(Enrichment_p < 0.1, na.rm=TRUE),
  n_significant_p001=sum(Enrichment_p < 0.01, na.rm=TRUE),
  median_enrichment=median(Enrichment, na.rm=TRUE),
  mean_enrichment=mean(Enrichment, na.rm=TRUE)
), by=.(view,dataset,celltype,trait,trait_label,method,method_label,panel,interpretation)]
fwrite(counts, file.path(paths$tables, "sldsc_count_table.tsv"), sep="\t")

paired <- x[matched_trait & is.finite(Enrichment) & Enrichment > 0,
  .(n=.N, mean_log2_enrichment=mean(log2(Enrichment)),
    median_log2_enrichment=median(log2(Enrichment)),
    mean_enrichment=mean(Enrichment), median_enrichment=median(Enrichment),
    significant_p01=sum(Enrichment_p < 0.1, na.rm=TRUE)),
  by=.(view,method,method_label)]
fwrite(paired, file.path(paths$tables, "sldsc_method_summary.tsv"), sep="\t")

theme_sldsc <- theme_bw(base_size=11) +
  theme(panel.grid.minor=element_blank(), axis.text.x=element_text(angle=35, hjust=1),
        strip.text=element_text(size=9), legend.position="bottom")

plot_trait <- function(z, view_name, trait_name) {
  z <- z[view == view_name & trait == trait_name]
  if (!nrow(z)) return(NULL)
  if (view_name == "celltype") {
    dodge <- position_dodge(width=.55)
    p <- ggplot(z, aes(method_label, Enrichment, color=celltype, group=celltype)) +
      geom_hline(yintercept=1, linetype=2, color="grey55") +
      geom_errorbar(aes(ymin=enrichment_low, ymax=enrichment_high), width=.12, position=dodge) +
      geom_point(position=dodge, size=2) +
      geom_text(aes(y=enrichment_high, label=significance), position=dodge,
                vjust=-.25, show.legend=FALSE, size=3) +
      facet_wrap(~dataset, ncol=3, scales="free_y") +
      scale_color_discrete(name="Cell type")
  } else {
    p <- ggplot(z, aes(method_label, Enrichment, color=method_label)) +
      geom_hline(yintercept=1, linetype=2, color="grey55") +
      geom_errorbar(aes(ymin=enrichment_low, ymax=enrichment_high), width=.12) +
      geom_point(size=2.2) +
      geom_text(aes(y=enrichment_high, label=significance), vjust=-.25,
                show.legend=FALSE, size=3) +
      facet_wrap(~dataset, ncol=3, scales="free_y") +
      scale_color_manual(values=method_colors, drop=FALSE, guide="none")
  }
  p + scale_x_discrete(drop=FALSE) +
    labs(title=paste0(z$trait_label[1], " — ", view_name, " annotations"),
         subtitle="S-LDSC: baseline-LD v2.2 + all tested peaks + raw-P<0.05 linked peaks; error bars ±1 SE; * P<0.1, ** P<0.01",
         x=NULL, y="Heritability enrichment") + theme_sldsc
}

matched <- x[matched_trait == TRUE]
for (v in c("pooled", "celltype")) {
  for (tr in unique(matched[view == v, trait])) {
    p <- plot_trait(matched, v, tr)
    if (!is.null(p)) {
      ggsave(file.path(paths$figures, paste0("sldsc_", tr, "_", v, ".pdf")),
             p, width=if (v == "celltype") 13 else 12, height=8.5)
    }
  }
}

# Compact method-level overview. This is secondary to the per-trait estimates;
# no comparator method is treated as a baseline.
p_summary <- ggplot(paired, aes(method_label, mean_log2_enrichment, color=method_label)) +
  geom_hline(yintercept=0, linetype=2, color="grey55") + geom_point(size=3) +
  geom_text(aes(label=paste0("n=", n)), vjust=-.8, size=3, show.legend=FALSE) +
  facet_wrap(~view) + scale_color_manual(values=method_colors, guide="none", drop=FALSE) +
  scale_x_discrete(drop=FALSE) +
  labs(title="Matched-trait S-LDSC summary", subtitle="Unweighted mean across available dataset × cell-type × trait estimates",
       x=NULL, y=expression(mean~log[2](heritability~enrichment))) + theme_sldsc
ggsave(file.path(paths$figures, "sldsc_matched_method_summary.pdf"), p_summary, width=10, height=5.5)

# Pituitary has no locally available matched GWAS; retain all six traits in a
# clearly labelled exploratory supplementary figure.
pit <- x[dataset == "Pituitary" & view == "pooled"]
if (nrow(pit)) {
  p_pit <- ggplot(pit, aes(method_label, Enrichment, color=method_label)) +
    geom_hline(yintercept=1, linetype=2, color="grey55") +
    geom_errorbar(aes(ymin=enrichment_low, ymax=enrichment_high), width=.12) + geom_point(size=2) +
    facet_wrap(~trait_label, scales="free_y", ncol=3) +
    scale_color_manual(values=method_colors, guide="none", drop=FALSE) + scale_x_discrete(drop=FALSE) +
    labs(title="Pituitary pooled annotations — exploratory traits",
         subtitle="No pituitary-matched GWAS summary statistic is present in the local reference panel",
         x=NULL, y="Heritability enrichment") + theme_sldsc
  ggsave(file.path(paths$figures, "sldsc_pituitary_exploratory_pooled.pdf"), p_pit, width=13, height=8)
}

message("Figures: ", paths$figures)
message("Tables: ", paths$tables)
