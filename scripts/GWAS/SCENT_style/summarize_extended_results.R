#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(if(length(args)) args[[1]] else ".",mustWork=TRUE)
base <- file.path(root,"results","validation","GWAS","SCENT_style")
sources <- c("FinnGen","RA","T1D","IBD")
x <- rbindlist(lapply(sources,function(s)
  fread(file.path(base,s,"tables","scent_summary.tsv"))),fill=TRUE)
fwrite(x,file.path(base,"combined_scent_summary.tsv.gz"),sep="\t")

immune_datasets <- c("Arthritis","pbmc","bmmc","Dogma_CTRL","Dogma_STIM","NEAT","SHARE")
metrics <- x[scope=="pooled" & dataset %chin% immune_datasets & is.finite(overall_enrichment),
  .(n_dataset_units=.N,mean_log2_1p_enrichment=mean(display_score),
    median_enrichment=median(overall_enrichment),mean_enrichment=mean(overall_enrichment)),
  by=.(source,method_label)]
setorder(metrics,source,-mean_log2_1p_enrichment)
fwrite(metrics,file.path(base,"pooled_immune_method_summary.tsv"),sep="\t")

wins <- x[scope=="pooled" & dataset %chin% immune_datasets & is.finite(overall_enrichment)]
wins[,is_best:=overall_enrichment==max(overall_enrichment),by=.(source,dataset)]
wins <- wins[is_best==TRUE,.(n_best=.N),by=.(source,method_label)]
setorder(wins,source,-n_best)
fwrite(wins,file.path(base,"pooled_immune_method_wins.tsv"),sep="\t")
