#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(if (length(args)) args[[1]] else ".", mustWork=TRUE)
base <- file.path(root,"results","validation","GWAS","SCENT_style")
source(file.path(root, "scripts", "GWAS", "SCENT_style", "source_results_paths.R"))
source_results <- source_results_paths(root)
restricted <- file.path(base,"testing","restricted_v2")
keep_methods <- c("epilink_poisson","epilink_nb","scmultimap",
                  "signac_corrected","signac_uncorrected")

read_scope <- function(scope) {
  f <- file.path(restricted,scope,"annotation_manifest.tsv")
  x <- fread(f)[method %chin% keep_methods]
  x[, `:=`(
    scope=scope,
    bed_path=file.path(restricted,scope,"SCENT","annotations",paste0(annotation_id,".hg38.bed")),
    correction=ifelse(grepl("^epilink",method),FALSE,
                      ifelse(method=="signac_corrected",TRUE,
                             ifelse(method=="signac_uncorrected",FALSE,NA)))
  )]
  x
}

ann <- rbindlist(list(read_scope("celltype"),read_scope("pooled")),fill=TRUE)

# Brain2 ALL/AD and Pituitary were prepared for the UKBB benchmark with the
# requested condition-effect EpiLink results and both Signac definitions.
b_pool <- fread(file.path(source_results$pooled_brain2,"manifest.tsv"))
b_pool <- b_pool[method %chin% keep_methods & dataset %chin% c("Brain2_ALL","Brain2_AD","Pituitary")]
b_pool[, `:=`(
  scope = "pooled",
  bed_path = file.path(source_results$pooled_annotations,paste0(annotation_id,".hg38.bed"))
)]
b_cell <- fread(file.path(source_results$celltype_brain2,"manifest.tsv"))
b_cell <- b_cell[method %chin% keep_methods & dataset %chin% c("Brain2_ALL","Brain2_AD","Pituitary")]
b_cell[, `:=`(
  scope = "celltype",
  bed_path = file.path(source_results$celltype_annotations,paste0(annotation_id,".hg38.bed"))
)]
b <- rbindlist(list(b_pool,b_cell),fill=TRUE)
ann <- rbindlist(list(ann,b),fill=TRUE,use.names=TRUE)

ann[, method_label := factor(method, levels=keep_methods,
  labels=c("EpiLink-Poisson","EpiLink-NB","scMultiMap",
           "Signac (corrected)","Signac (uncorrected)"))]
ann[, method_label := as.character(method_label)]
ann <- ann[file.exists(bed_path)]
ann <- unique(ann,by=c("annotation_id","scope"))
setorder(ann,scope,dataset,celltype,method)

out <- file.path(base,"config","extended_annotation_manifest.tsv")
fwrite(ann,out,sep="\t")
colors <- data.table(
  method=keep_methods,
  method_label=c("EpiLink-Poisson","EpiLink-NB","scMultiMap",
                 "Signac (corrected)","Signac (uncorrected)"),
  color=c("#E47159","#C87878","#edd089","#54b2a2","#78A878")
)
fwrite(colors,file.path(base,"config","extended_method_colors.tsv"),sep="\t")
message("Annotations: ",nrow(ann)," (celltype=",sum(ann$scope=="celltype"),
        ", pooled=",sum(ann$scope=="pooled"),")")
message("Datasets: ",paste(unique(ann$dataset),collapse=", "))
