#!/usr/bin/env Rscript

# Screen Brain2 AD-vs-control differential EpiLink-Poisson pairs for overlap
# with fine-mapped neurological GWAS variants.  scMultiMap-diff is the direct
# differential comparator.  Signac has no differential output in this project,
# so its ALL/AD/CTRL calls are reported only as an auxiliary exclusion column.

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[[1L]] else ".", mustWork = TRUE)
out_dir <- file.path(root, "results", "validation", "functional", "GWAS_loci")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_two <- function(path) {
  if (!file.exists(path)) return(data.table(gene = character(), peak = character()))
  unique(fread(path, select = c("gene", "peak"), showProgress = FALSE))
}

read_gwas <- function() {
  ref <- file.path(root, "results", "validation", "GWAS", "reference")
  out <- list()
  for (src in c("FinnGen", "RA", "IBD", "T1D")) {
    p <- file.path(ref, src, "processed", "causal_variants_PIPgt0.2.tsv.gz")
    if (file.exists(p)) out[[src]] <- fread(p, showProgress = FALSE)
  }
  uk <- file.path(ref, "UKBB", "ukbb_94traits_causal.hg38.tsv.gz")
  if (file.exists(uk)) {
    z <- fread(uk, showProgress = FALSE)
    z[, `:=`(source = "UKBB", pip = NA_real_,
             locus = paste(trait, chr, pos, sep = "::"),
             consequence = NA_character_,
             variant_id = paste0(sub("^chr", "", chr), ":", pos))]
    out$UKBB <- z
  }
  x <- rbindlist(out, fill = TRUE, use.names = TRUE)
  x[, variant_chr := ifelse(grepl("^chr", chr), chr, paste0("chr", chr))]
  x[, chr := NULL]
  rx <- "alzheimer|dement|cognit|parkinson|schizo|bipolar|depress|neurotic|insomnia|migraine"
  unique(x[grepl(rx, trait, ignore.case = TRUE)],
         by = c("source", "trait", "variant_chr", "pos", "variant_id"))
}

configs <- data.table(
  celltype = c("Astrocytes", "Microglia"),
  ep_dir = c("epilink_diff_poisson", "epilink_diff_poisson_donor_effect_donor_robust")
)

all_hits <- list()
for (i in seq_len(nrow(configs))) {
  ct <- configs$celltype[i]
  base <- file.path(root, "results", "real", "Brain2_diff", ct)
  ep_path <- file.path(base, configs$ep_dir[i], "all_pairs.tsv.gz")
  if (!file.exists(ep_path)) next
  ep <- fread(ep_path, select = c("gene", "peak", "beta_diff", "pval", "padj",
                                  "calibrated_q", "significant"), showProgress = FALSE)
  ep <- ep[is.finite(padj) & padj < 0.1]
  ep[, `:=`(celltype = ct, epilink_model = configs$ep_dir[i])]

  sm <- read_two(file.path(base, "scmultimap", "significant_pairs.tsv.gz"))
  sm[, scmultimap_diff := TRUE]
  ep <- merge(ep, sm, by = c("gene", "peak"), all.x = TRUE)
  ep[is.na(scmultimap_diff), scmultimap_diff := FALSE]

  sg <- list()
  for (condition in c("ALL", "AD", "CTRL")) {
    sdir <- file.path(root, "results", "real", "Brain2", ct, condition, "signac")
    for (f in c("significant_pairs.tsv.gz", "significant_pairs_uncorrected.tsv.gz"))
      sg[[length(sg) + 1L]] <- read_two(file.path(sdir, f))
  }
  sg <- unique(rbindlist(sg, fill = TRUE))
  sg[, signac_any_state := TRUE]
  ep <- merge(ep, sg, by = c("gene", "peak"), all.x = TRUE)
  ep[is.na(signac_any_state), signac_any_state := FALSE]

  bits <- tstrsplit(ep$peak, "-", fixed = TRUE)
  ep[, `:=`(peak_chr = bits[[1L]], peak_start = as.integer(bits[[2L]]),
            peak_end = as.integer(bits[[3L]]))]
  setkey(ep, peak_chr, peak_start, peak_end)
  gv <- copy(read_gwas())
  gv[, `:=`(variant_start = as.integer(pos), variant_end = as.integer(pos))]
  setkey(gv, variant_chr, variant_start, variant_end)
  hit <- foverlaps(gv, ep,
                   by.x = c("variant_chr", "variant_start", "variant_end"),
                   by.y = c("peak_chr", "peak_start", "peak_end"),
                   type = "within", nomatch = 0L)
  if (nrow(hit)) all_hits[[length(all_hits) + 1L]] <- hit
}

ans <- rbindlist(all_hits, fill = TRUE)
if (nrow(ans)) {
  ans[, strict_diff_epilink_only := !scmultimap_diff]
  ans[, auxiliary_signac_absent := !signac_any_state]
  setorder(ans, -strict_diff_epilink_only, -auxiliary_signac_absent, -pip, padj)
}
fwrite(ans, file.path(out_dir, "brain2_diff_poisson_gwas_candidates.tsv.gz"), sep = "\t")
message("Wrote ", nrow(ans), " Brain2 differential GWAS-overlap rows")
