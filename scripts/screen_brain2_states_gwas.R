#!/usr/bin/env Rscript

# Strict state-matched Brain2 screen: EpiLink-Poisson FDR < 0.1 and the exact
# peak-gene pair absent from scMultiMap, corrected Signac and uncorrected Signac
# calls in the same cell type and condition (ALL, AD or CTRL).

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

gwas <- read_gwas()
hits <- list()
for (ct in c("Astrocytes", "Microglia")) for (state in c("ALL", "AD", "CTRL")) {
  base <- file.path(root, "results", "real", "Brain2", ct, state)
  p <- file.path(base, "epilink_poisson", "all_pairs.tsv.gz")
  if (!file.exists(p) && state == "ALL")
    p <- file.path(base, "epilink_poisson_condition_effect", "all_pairs.tsv.gz")
  if (!file.exists(p)) next
  hdr <- names(fread(p, nrows = 0L, showProgress = FALSE))
  cols <- intersect(c("gene", "peak", "beta", "pval", "pvalue", "padj"), hdr)
  ep <- fread(p, select = cols, showProgress = FALSE)
  ep <- ep[is.finite(padj) & padj < 0.1]
  if (!nrow(ep)) next
  ep[, `:=`(celltype = ct, state = state)]

  sm <- read_two(file.path(base, "scmultimap", "significant_pairs.tsv.gz"))
  sm[, scmultimap := TRUE]
  sg <- rbindlist(list(
    read_two(file.path(base, "signac", "significant_pairs.tsv.gz")),
    read_two(file.path(base, "signac", "significant_pairs_uncorrected.tsv.gz"))
  ))
  sg <- unique(sg)
  sg[, signac := TRUE]
  ep <- merge(ep, sm, by = c("gene", "peak"), all.x = TRUE)
  ep <- merge(ep, sg, by = c("gene", "peak"), all.x = TRUE)
  ep[is.na(scmultimap), scmultimap := FALSE]
  ep[is.na(signac), signac := FALSE]
  ep[, epilink_only := !scmultimap & !signac]

  bits <- tstrsplit(ep$peak, "-", fixed = TRUE)
  ep[, `:=`(peak_chr = bits[[1L]], peak_start = as.integer(bits[[2L]]),
            peak_end = as.integer(bits[[3L]]))]
  setkey(ep, peak_chr, peak_start, peak_end)
  gv <- copy(gwas)
  gv[, `:=`(variant_start = as.integer(pos), variant_end = as.integer(pos))]
  setkey(gv, variant_chr, variant_start, variant_end)
  z <- foverlaps(gv, ep,
                 by.x = c("variant_chr", "variant_start", "variant_end"),
                 by.y = c("peak_chr", "peak_start", "peak_end"),
                 type = "within", nomatch = 0L)
  if (nrow(z)) hits[[length(hits) + 1L]] <- z
}

ans <- rbindlist(hits, fill = TRUE)
if (nrow(ans)) setorder(ans, -epilink_only, -pip, padj, celltype, state, gene)
fwrite(ans, file.path(out_dir, "brain2_state_matched_poisson_gwas_candidates.tsv.gz"), sep = "\t")
message("Wrote ", nrow(ans), " state-matched rows; ",
        if (nrow(ans)) ans[, sum(epilink_only)] else 0L, " EpiLink-only")
