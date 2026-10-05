#!/usr/bin/env Rscript

# Check Yang Supplementary Table 6 matches among non-differential Brain2
# Microglia significant peak-gene pairs.

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
project_dir <- if (length(args)) normalizePath(args[[1L]]) else normalizePath(".")
setwd(project_dir)

match_file <- file.path(project_dir, "results", "validation", "functional",
                        "Brain2_diff", "Microglia", "yang_overlap",
                        "yang_brain2_peak_gene_matches.tsv.gz")
if (!file.exists(match_file)) {
  stop("Run scripts/analyze_brain2_microglia_yang.R first: missing ", match_file,
       call. = FALSE)
}

out_dir <- file.path(project_dir, "results", "validation", "functional",
                     "Brain2", "Microglia", "yang_overlap_nondiff")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

matches <- fread(match_file, showProgress = FALSE)
matches[, `:=`(gene = as.character(gene), peak = as.character(peak))]

methods <- c(
  epilink_poisson = file.path(project_dir, "results", "real", "Brain2", "Microglia", "epilink_poisson"),
  epilink_nb = file.path(project_dir, "results", "real", "Brain2", "Microglia", "epilink_nb"),
  poisson_glm = file.path(project_dir, "results", "real", "Brain2", "Microglia", "poisson_glm"),
  nb_glm = file.path(project_dir, "results", "real", "Brain2", "Microglia", "nb_glm"),
  scmultimap = file.path(project_dir, "results", "real", "Brain2", "Microglia", "scmultimap")
)

summary <- rbindlist(lapply(names(methods), function(method) {
  sig_file <- file.path(methods[[method]], "significant_pairs.tsv.gz")
  sig <- fread(sig_file, showProgress = FALSE)
  if (!all(c("gene", "peak") %in% names(sig))) {
    sig <- data.table(gene = character(), peak = character())
  }
  sig[, `:=`(gene = as.character(gene), peak = as.character(peak))]
  hit <- merge(matches, sig, by = c("gene", "peak"), all = FALSE,
               allow.cartesian = TRUE)
  fwrite(hit, file.path(out_dir, paste0(method, "_yang_significant_matches.tsv.gz")),
         sep = "\t")
  data.table(
    method = method,
    significant_pairs = uniqueN(sig[, .(gene, peak)]),
    exact_yang_matches = uniqueN(hit[, .(gene, peak)]),
    Yang_variant_target_rows = nrow(hit)
  )
}))

fwrite(summary, file.path(out_dir, "summary.tsv"), sep = "\t")
message("Wrote non-differential Yang matching results to ", normalizePath(out_dir))
print(summary)
