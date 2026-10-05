#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
root <- "results/validation/eQTL/GTEx/eQTL1"
summaries <- Sys.glob(file.path(root, "*", "*", "**", "summary.tsv"))
# Sys.glob does not implement recursive ** consistently; use list.files and
# retain method-level summaries only.
summaries <- list.files(root, pattern = "^summary\\.tsv$", recursive = TRUE,
                        full.names = TRUE)
summaries <- summaries[grepl("/(epilink_poisson|epilink_nb|scmultimap|signac_corrected|signac_uncorrected)/summary\\.tsv$",
                             summaries)]
if (!length(summaries)) stop("No eQTL1 method summaries found")
x <- rbindlist(lapply(summaries, fread), fill = TRUE)
x[, pooled_or_padj := p.adjust(pooled_or_pval, method = "BH")]
setorder(x, scope, dataset, celltype, method)
fwrite(x, file.path(root, "all_results.tsv"), sep = "\t")

tissue <- list.files(root, pattern = "^tissue_summary\\.tsv\\.gz$", recursive = TRUE,
                     full.names = TRUE)
t <- rbindlist(lapply(tissue, fread), fill = TRUE)
setorder(t, scope, dataset, celltype, method, tissue)
fwrite(t, file.path(root, "all_tissue_results.tsv.gz"), sep = "\t")

audits <- list.files(root, pattern = "^scope_audit\\.tsv$", recursive = TRUE,
                     full.names = TRUE)
a <- unique(rbindlist(lapply(audits, fread), fill = TRUE))
setorder(a, scope, dataset, celltype)
fwrite(a, file.path(root, "scope_audit.tsv"), sep = "\t")
message("Collected ", nrow(x), " primary summaries, ", nrow(t),
        " tissue summaries and ", nrow(a), " scope audits")
