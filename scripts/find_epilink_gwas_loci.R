#!/usr/bin/env Rscript

# Find fine-mapped GWAS variants in EpiLink links that are not significant in
# scMultiMap or Signac.  This is a candidate-generation step for locus-level
# functional follow-up, not a formal enrichment test.

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[[1L]] else ".", mustWork = TRUE)
out_dir <- file.path(root, "results", "validation", "functional", "GWAS_loci")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_pairs <- function(path, method, dataset, celltype, analysis = "standard") {
  if (!file.exists(path) || file.info(path)$size == 0) return(NULL)
  # all_pairs files can be large; read only columns needed for the screen.
  hdr <- tryCatch(names(fread(path, nrows = 0L, showProgress = FALSE)), error = function(e) character())
  keep <- intersect(c("gene", "peak", "padj", "pvalue", "pval", "score", "beta", "cor"), hdr)
  x <- tryCatch(fread(path, select = keep, showProgress = FALSE), error = function(e) NULL)
  if (is.null(x) || !nrow(x) || !all(c("gene", "peak") %chin% names(x))) return(NULL)
  stat <- intersect(c("padj", "pvalue", "pval", "score", "beta", "cor"), names(x))
  x <- x[, c("gene", "peak", stat), with = FALSE]
  if (method == "epilink_poisson" && "padj" %chin% names(x))
    x <- x[is.finite(padj) & padj < 0.1]
  x[, `:=`(method = method, dataset = dataset, celltype = celltype,
           analysis = analysis)]
  unique(x, by = c("gene", "peak", "method", "dataset", "celltype", "analysis"))
}

result_path <- function(dataset, celltype, method) {
  base <- file.path(root, "results", "real", dataset, celltype)
  if (dataset == "Brain2") {
    base <- file.path(base, "ALL")
    choices <- switch(method,
      epilink_poisson = c("epilink_poisson_condition_effect", "epilink_poisson", "epilink_poisson_donor_effect"),
      scmultimap = "scmultimap", signac = "signac")
    dirs <- file.path(base, choices)
    dirs <- dirs[file.exists(dirs)]
    if (!length(dirs)) return(NA_character_)
    base <- dirs[[1L]]
  } else {
    base <- file.path(base, method)
  }
  if (method == "epilink_poisson") file.path(base, "all_pairs.tsv.gz") else
    file.path(base, "significant_pairs.tsv.gz")
}

panel <- list(
  immune = list(
    datasets = c("Arthritis", "Dogma_CTRL", "Dogma_STIM", "NEAT", "SHARE", "bmmc", "pbmc"),
    trait_regex = paste(c("rheuma", "arthritis", "crohn", "colitis", "bowel", "IBD",
                          "T1D", "type 1 diabetes", "multiple.sclerosis", "AUTOIMMUNE",
                          "asthma", "atopic", "hypothy"), collapse = "|")),
  pituitary = list(
    datasets = "Pituitary",
    trait_regex = paste(c("menarch", "menop", "endometri", "menorr", "ovarian", "ovary",
                          "uter", "hypogonad", "fertil", "FSH", "LH", "testosterone",
                          "SHBG", "IGF1", "height"), collapse = "|")),
  brain = list(
    datasets = c("Brain", "Brain2"),
    trait_regex = "alzheimer|dement|cognit|parkinson|schizo|bipolar|depress|neurotic|insomnia|migraine")
)

all_links <- list(); k <- 0L
for (tissue in names(panel)) {
  for (dataset in panel[[tissue]]$datasets) {
    droot <- file.path(root, "results", "real", dataset)
    if (!dir.exists(droot)) next
    celltypes <- list.dirs(droot, recursive = FALSE, full.names = FALSE)
    for (celltype in celltypes) {
      for (method in c("epilink_poisson", "scmultimap", "signac")) {
        path <- result_path(dataset, celltype, method)
        if (is.na(path)) next
        z <- read_pairs(path, method, dataset, celltype,
                        if (dataset == "Brain2") basename(dirname(path)) else "standard")
        if (!is.null(z)) {
          k <- k + 1L
          z[, tissue := tissue]
          all_links[[k]] <- z
        }
      }
    }
  }
}
links <- rbindlist(all_links, fill = TRUE)
if (!nrow(links)) stop("No significant result pairs found")

parse_peak <- function(x) {
  bits <- tstrsplit(x, "-", fixed = TRUE)
  data.table(peak_chr = bits[[1L]], peak_start = as.integer(bits[[2L]]),
             peak_end = as.integer(bits[[3L]]))
}
links <- cbind(links, parse_peak(links$peak))
links[, pair_key := paste(dataset, celltype, gene, peak, sep = "|")]

# Exact-pair competitor status. Signac's exported significant table is the
# corrected primary call; the uncorrected table is additionally checked below.
competitors <- unique(links[method %chin% c("scmultimap", "signac"),
                            .(pair_key, competitor = method)])
wide <- dcast(competitors, pair_key ~ competitor, fun.aggregate = length,
              value.var = "competitor")
for (nm in c("scmultimap", "signac")) if (!nm %chin% names(wide)) wide[, (nm) := 0L]

# Include Signac's uncorrected calls, making the EpiLink-only criterion strict.
uncorrected <- list(); u <- 0L
for (tissue in names(panel)) for (dataset in panel[[tissue]]$datasets) {
  droot <- file.path(root, "results", "real", dataset)
  if (!dir.exists(droot)) next
  for (celltype in list.dirs(droot, recursive = FALSE, full.names = FALSE)) {
    sdir <- if (dataset == "Brain2") file.path(droot, celltype, "ALL", "signac") else
      file.path(droot, celltype, "signac")
    p <- file.path(sdir, "significant_pairs_uncorrected.tsv.gz")
    if (!file.exists(p)) next
    x <- tryCatch(fread(p, select = c("gene", "peak"), showProgress = FALSE),
                  error = function(e) NULL)
    if (!is.null(x) && nrow(x)) {
      u <- u + 1L
      x[, pair_key := paste(dataset, celltype, gene, peak, sep = "|")]
      uncorrected[[u]] <- unique(x[, .(pair_key)])
    }
  }
}
uncorrected_keys <- if (length(uncorrected)) unique(rbindlist(uncorrected)$pair_key) else character()

epi <- links[method == "epilink_poisson" & is.finite(padj) & padj < 0.1]
epi <- merge(epi, wide, by = "pair_key", all.x = TRUE)
epi[is.na(scmultimap), scmultimap := 0L]
epi[is.na(signac), signac := 0L]
epi[, signac_uncorrected := pair_key %chin% uncorrected_keys]
epi[, epilink_only := scmultimap == 0L & signac == 0L & !signac_uncorrected]

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
    z[, `:=`(source = "UKBB", pip = NA_real_, locus = paste(trait, chr, pos, sep = "::"),
             consequence = NA_character_, variant_id = paste0(sub("^chr", "", chr), ":", pos))]
    out$UKBB <- z
  }
  ans <- rbindlist(out, fill = TRUE, use.names = TRUE)
  ans[, variant_chr := ifelse(grepl("^chr", chr), chr, paste0("chr", chr))]
  ans[, chr := NULL]
  unique(ans, by = c("source", "trait", "variant_chr", "pos", "variant_id"))
}
gwas <- read_gwas()

# Exact interval overlap between fine-mapped variant positions and linked peaks.
hits <- list(); h <- 0L
for (tissue_name in names(panel)) {
  e <- epi[tissue == tissue_name]
  v <- gwas[grepl(panel[[tissue_name]]$trait_regex, trait, ignore.case = TRUE)]
  if (!nrow(e) || !nrow(v)) next
  ee <- copy(e)
  vv <- copy(v)
  setkey(ee, peak_chr, peak_start, peak_end)
  vv[, `:=`(variant_start = as.integer(pos), variant_end = as.integer(pos))]
  setkey(vv, variant_chr, variant_start, variant_end)
  ov <- foverlaps(vv, ee,
                  by.x = c("variant_chr", "variant_start", "variant_end"),
                  by.y = c("peak_chr", "peak_start", "peak_end"),
                  type = "within", nomatch = 0L)
  if (nrow(ov)) {
    h <- h + 1L
    hits[[h]] <- ov
  }
}
candidates <- if (length(hits)) rbindlist(hits, fill = TRUE) else data.table()
if (nrow(candidates)) {
  candidates[, peak_mid := (peak_start + peak_end) / 2]
  candidates[, variant_peak_offset := pos - peak_mid]
  setorder(candidates, tissue, -epilink_only, -pip, dataset, celltype, gene)
}

fwrite(candidates, file.path(out_dir, "epilink_finemapped_candidates.tsv.gz"), sep = "\t")
fwrite(candidates[epilink_only == TRUE],
       file.path(out_dir, "epilink_only_finemapped_candidates.tsv.gz"), sep = "\t")
summary <- candidates[, .(
  n_variant_peak_gene = .N,
  n_variants = uniqueN(paste(source, trait, variant_id)),
  n_genes = uniqueN(gene),
  n_epilink_only = sum(epilink_only)
), by = .(tissue, dataset, celltype, method)]
fwrite(summary, file.path(out_dir, "candidate_summary.tsv"), sep = "\t")

manifest <- data.table(
  item = c("fine_mapping", "EpiLink call", "competitor exclusion", "genome build"),
  definition = c("PIP > 0.2 for FinnGen/RA/IBD/T1D; published PIP-filtered UKBB causal table",
                 "EpiLink-Poisson all_pairs.tsv.gz filtered at ordinary FDR < 0.1; NB excluded",
                 "exact peak-gene pair absent from scMultiMap, corrected Signac, and uncorrected Signac significant calls in the same dataset/cell type",
                 "GRCh38")
)
fwrite(manifest, file.path(out_dir, "candidate_manifest.tsv"), sep = "\t")
message("Wrote ", nrow(candidates), " candidate rows (",
        candidates[, sum(epilink_only)], " EpiLink-only) to ", out_dir)
