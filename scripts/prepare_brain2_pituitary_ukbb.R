#!/usr/bin/env Rscript

# Build UKBB/SCENT-style annotations for Brain2 ALL/CTRL/AD and Pituitary.
# This is separate from the original manifest because Brain2 has state-specific
# result directories and ALL Microglia uses donor-effect results.

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
})

args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[[1]] else ".", mustWork = TRUE)
fdr <- if (length(args) >= 2L) as.numeric(args[[2]]) else 0.1
gwas <- file.path(root, "results", "validation", "GWAS")
base <- file.path(gwas, "SCENT", "sakaue_style", "brain2_pituitary_ukbb")
ann_dir <- file.path(base, "annotations")
dir.create(ann_dir, recursive = TRUE, showWarnings = FALSE)
source(file.path(root, "R", "EpiLink.R"))

safe_id <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

read_peak_table <- function(path) {
  if (!file.exists(path)) return(data.table(peak = character()))
  tryCatch(fread(path, select = "peak", showProgress = FALSE),
           error = function(e) data.table(peak = character()))
}

parse_peaks <- function(x) {
  x <- unique(x[!is.na(x) & nzchar(x)])
  if (!length(x)) return(data.table(chr = character(), start = integer(), end = integer()))
  z <- tstrsplit(x, "-", fixed = TRUE)
  if (length(z) != 3L) stop("Peak IDs must have chr-start-end form")
  ans <- data.table(chr = z[[1]], start = as.integer(z[[2]]), end = as.integer(z[[3]]))
  ans <- ans[grepl("^chr([1-9]|1[0-9]|2[0-2]|X|Y)$", chr) &
               !is.na(start) & !is.na(end) & start <= end]
  unique(ans)
}

reduce_bed <- function(x) {
  if (!nrow(x)) return(data.table(chr = character(), start = integer(), end = integer()))
  gr <- reduce(GRanges(x$chr, IRanges(x$start, x$end)), ignore.strand = TRUE)
  data.table(chr = as.character(seqnames(gr)), start = start(gr), end = end(gr))
}

write_bed <- function(x, path) {
  y <- copy(x)
  y[, start0 := pmax(start - 1L, 0L)]
  fwrite(y[, .(chr, start0, end)], path, sep = "\t", col.names = FALSE)
}

scope <- data.table(
  dataset = c("Brain2_ALL", "Brain2_CTRL", "Brain2_AD", "Pituitary"),
  source_dataset = c("Brain2", "Brain2", "Brain2", "Pituitary"),
  state = c("ALL", "CTRL", "AD", ""),
  panel = c("brain", "brain", "brain", "pituitary"),
  celltypes = c("Astrocytes;Microglia", "Astrocytes;Microglia",
                "Astrocytes;Microglia",
                "Corticotropes;Gonadotropes;Lactotropes;Somatotropes;Stem cells;Thyrotropes")
)
scope[, celltype_vec := strsplit(celltypes, ";", fixed = TRUE)]

methods <- data.table(
  method = c("epilink_nb", "epilink_poisson", "scmultimap",
             "signac_corrected", "signac_uncorrected"),
  method_label = c("EpiLink-NB", "EpiLink-Poisson", "scMultiMap",
                   "Signac (corrected)", "Signac (uncorrected)"),
  correction = c(NA, NA, NA, TRUE, FALSE)
)

# Prefer ordinary association results. For Brain2 ALL Microglia, ordinary
# EpiLink was not run, so fall back to donor-effect association results.
method_dir <- function(source_dataset, state, celltype, method) {
  root_dir <- file.path(root, "results", "real", source_dataset, celltype)
  if (source_dataset == "Brain2") {
    candidates <- if (state == "ALL") {
      if (method == "epilink_nb") c("epilink_nb_condition_effect", "epilink_nb", "epilink_nb_donor_effect")
      else if (method == "epilink_poisson") c("epilink_poisson_condition_effect", "epilink_poisson", "epilink_poisson_donor_effect")
      else if (method == "scmultimap") c("scmultimap", "scmultimap_donor_effect")
      else c("signac")
    } else {
      if (method == "epilink_nb") c("epilink_nb_condition_effect", "epilink_nb")
      else if (method == "epilink_poisson") c("epilink_poisson_condition_effect", "epilink_poisson", "epilink_poisson_donor_effect")
      else if (method == "scmultimap") c("scmultimap", "scmultimap_donor_effect")
      else c("signac")
    }
    dirs <- file.path(root_dir, state, candidates)
  } else {
    candidates <- if (method %chin% c("epilink_nb", "epilink_poisson")) method else "signac"
    dirs <- file.path(root_dir, candidates)
  }
  dirs[file.exists(dirs)]
}

extract_peaks <- function(source_dataset, state, celltype, method) {
  dirs <- method_dir(source_dataset, state, celltype, method)
  if (!length(dirs)) return(list(peaks = character(), result_dir = NA_character_,
                                 variant = "missing", correction = NA))
  d <- dirs[[1L]]
  if (method %chin% c("epilink_nb", "epilink_poisson")) {
    # Use the completed exported table, avoiding repeated loads of large RDS objects.
    p <- file.path(d, "significant_pairs.tsv.gz")
    md <- fread(file.path(d, "run_metadata.tsv"), nrows = 1L, showProgress = FALSE)
    actual_correction <- if ("correction" %in% names(md)) as.logical(md$correction[1L]) else NA
    variant <- if ("condition_effect" %in% names(md) && isTRUE(as.logical(md$condition_effect[1L])))
      "condition_effect" else if ("donor_effect" %in% names(md) && isTRUE(as.logical(md$donor_effect[1L])))
      "donor_effect" else if ("condition_subset" %in% names(md) && !is.na(md$condition_subset[1L]))
      "condition_subset" else "standard"
    return(list(peaks = read_peak_table(p)$peak, result_dir = d,
                variant = variant, correction = actual_correction))
  }
  if (method == "signac_corrected") {
    return(list(peaks = read_peak_table(file.path(d, "significant_pairs_corrected.tsv.gz"))$peak,
                result_dir = d, variant = "corrected", correction = TRUE))
  }
  if (method == "signac_uncorrected") {
    return(list(peaks = read_peak_table(file.path(d, "significant_pairs_uncorrected.tsv.gz"))$peak,
                result_dir = d, variant = "uncorrected", correction = FALSE))
  }
  list(peaks = read_peak_table(file.path(d, "significant_pairs.tsv.gz"))$peak,
       result_dir = d, variant = "reported_significant_pairs", correction = NA)
}

rows <- list(); k <- 0L
for (i in seq_len(nrow(scope))) {
  sc <- scope[i]
  for (ct in sc$celltype_vec[[1L]]) {
    for (j in seq_len(nrow(methods))) {
      mi <- methods[j]
      got <- extract_peaks(sc$source_dataset, sc$state, ct, mi$method)
      native <- reduce_bed(parse_peaks(got$peaks))
      id <- paste(safe_id(sc$dataset), safe_id(ct), safe_id(mi$method), sep = "__")
      write_bed(native, file.path(ann_dir, paste0(id, ".hg38.bed")))
      k <- k + 1L
      rows[[k]] <- data.table(
        annotation_id = id, dataset = sc$dataset, celltype = ct,
        celltype_scope = "celltype_specific", n_celltypes = 1L, celltypes = ct,
        panel = sc$panel, interpretation = "exploratory", method = mi$method,
        method_label = mi$method_label,
        correction = ifelse(is.na(got$correction), mi$correction, got$correction), fdr = fdr,
        source_build = "hg38", n_pairs = NA_integer_,
        n_native_intervals = nrow(native), n_hg38_intervals = nrow(native),
        n_hg19_intervals = NA_integer_, lift_fraction_hg38 = 1,
        lift_fraction_hg19 = NA_real_, result_dir = got$result_dir,
        analysis_variant = got$variant
      )
    }
  }
  for (j in seq_len(nrow(methods))) {
    mi <- methods[j]
    extracted <- lapply(sc$celltype_vec[[1L]], function(ct)
      extract_peaks(sc$source_dataset, sc$state, ct, mi$method))
    pieces <- lapply(extracted, function(got) parse_peaks(got$peaks))
    native <- reduce_bed(rbindlist(pieces, fill = TRUE))
    id <- paste(safe_id(sc$dataset), "ALL_REQUESTED", safe_id(mi$method), sep = "__")
    write_bed(native, file.path(ann_dir, paste0(id, ".hg38.bed")))
    variants <- unique(vapply(extracted, function(got) got$variant, character(1)))
    corrections <- unique(vapply(extracted, function(got) {
      if (is.na(got$correction)) mi$correction else got$correction
    }, logical(1)))
    k <- k + 1L
    rows[[k]] <- data.table(
      annotation_id = id, dataset = sc$dataset, celltype = "ALL_REQUESTED",
      celltype_scope = "pooled", n_celltypes = length(sc$celltype_vec[[1L]]),
      celltypes = sc$celltypes, panel = sc$panel, interpretation = "exploratory",
      method = mi$method, method_label = mi$method_label,
      correction = if (length(corrections) == 1L) corrections else NA,
      fdr = fdr, source_build = "hg38", n_pairs = NA_integer_,
      n_native_intervals = nrow(native), n_hg38_intervals = nrow(native),
      n_hg19_intervals = NA_integer_, lift_fraction_hg38 = 1,
      lift_fraction_hg19 = NA_real_, result_dir = paste(variants, collapse = ";"),
      analysis_variant = "union"
    )
  }
}

manifest <- rbindlist(rows, fill = TRUE)
fwrite(manifest, file.path(base, "sakaue_annotation_manifest.tsv"), sep = "\t")
fwrite(manifest[, .(dataset, celltype, celltype_scope, method, method_label,
                    correction, n_native_intervals, n_hg38_intervals,
                    result_dir, analysis_variant)],
       file.path(base, "ukbb_annotation_counts.tsv"), sep = "\t")
fwrite(scope[, .(dataset, source_dataset, state, panel, celltypes)],
       file.path(base, "dataset_scope.tsv"), sep = "\t")
message("Wrote ", nrow(manifest), " annotations to ", base)
