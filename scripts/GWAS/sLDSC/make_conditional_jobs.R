#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(if (length(args)) args[[1]] else ".", mustWork=TRUE)
view <- if (length(args) >= 2L) args[[2]] else "celltype"
source(file.path(root, "scripts", "GWAS", "sLDSC", "source_results_paths.R"))
paths <- sldsc_source_results_paths(root)
base <- paths[[view]]
m <- fread(file.path(base, "manifest.tsv"))
t <- fread(file.path(root, "results", "validation", "GWAS", "SCENT_style", "config", "traits.tsv"))
t <- unique(t[, .(trait, trait_label=label)])

trait_set <- function(panel) {
  if (panel == "immune") return(c("M13_RHEUMA", "T1D", "K11_IBD_STRICT", "G6_MS"))
  if (panel == "arthritis_context") return("M13_RHEUMA")
  if (panel == "brain") return("G6_ALZHEIMER")
  # No disease-matched pituitary summary statistic is present locally. Keep
  # all six available traits as explicitly exploratory results.
  t$trait
}

jobs <- rbindlist(lapply(seq_len(nrow(m)), function(i) {
  tr <- trait_set(m$panel[i])
  data.table(annotation_id=m$annotation_id[i], background_id=m$background_id[i],
             panel=m$panel[i], trait=tr)
}))
fwrite(jobs, file.path(base, "ldscores", "sldsc_jobs.tsv"), sep="\t")
message(view, ": ", nrow(jobs), " conditional S-LDSC jobs")
