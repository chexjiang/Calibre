# Shared validation helpers for CRISPR, chromatin-contact, and genetic benchmarks.

# New Signac results written by scripts/run_real.R contain both the
# corrected and uncorrected p-values. Older results only contain the corrected
# pvalue column, so retain the /2 fallback for backward compatibility.
validation_signac_method_ids <- c("signac", "signac_uncorrected")

validation_result_method_id <- function(method_id) {
  if (identical(method_id, "signac_uncorrected")) "signac" else method_id
}

validation_signac_pvalues <- function(result, method_id) {
  if (!method_id %in% validation_signac_method_ids) {
    stop("Unsupported Signac validation method: ", method_id, call. = FALSE)
  }
  if (identical(method_id, "signac_uncorrected") &&
      "pvalue_uncorrected" %in% names(result)) {
    pvalue <- as.numeric(result$pvalue_uncorrected)
  } else if (identical(method_id, "signac") &&
             "pvalue_corrected" %in% names(result)) {
    pvalue <- as.numeric(result$pvalue_corrected)
  } else {
    if (!"pvalue" %in% names(result)) {
      stop("Signac result lacks pvalue", call. = FALSE)
    }
    pvalue <- as.numeric(result$pvalue)
    if (identical(method_id, "signac_uncorrected")) pvalue <- pvalue / 2
  }
  pmin(pmax(pvalue, 0), 1)
}

links_to_ranges <- function(links, coordinate_chain = NULL) {
  links <- unique(links[, c("gene", "peak"), drop = FALSE])
  peaks <- as.character(links$peak)
  parsed <- regexec("^(.*?)[-:]([0-9]+)[-:]([0-9]+)$", peaks)
  fields <- regmatches(peaks, parsed)
  valid <- lengths(fields) == 4L
  if (!all(valid)) {
    stop("unable to parse peak coordinates: ", paste(head(peaks[!valid], 3L), collapse = ", "), call. = FALSE)
  }
  ranges <- GenomicRanges::GRanges(
    seqnames = vapply(fields, `[[`, character(1), 2L),
    ranges = IRanges::IRanges(
      start = as.integer(vapply(fields, `[[`, character(1), 3L)),
      end = as.integer(vapply(fields, `[[`, character(1), 4L))
    )
  )
  if (!is.null(coordinate_chain)) {
    if (is.character(coordinate_chain)) {
      coordinate_chain <- rtracklayer::import.chain(coordinate_chain)
    }
    lifted <- rtracklayer::liftOver(ranges, coordinate_chain)
    keep <- lengths(lifted) == 1L
    ranges <- unlist(lifted[keep], use.names = FALSE)
    links <- links[keep, , drop = FALSE]
  }
  ranges$gene <- as.character(links$gene)
  ranges$pair_key <- paste(links$gene, links$peak, sep = "\r")
  ranges
}

read_interval_gene_validation <- function(path) {
  table <- utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("chr", "start", "end", "gene")
  if (!all(required %in% names(table))) stop("validation file lacks: ", paste(setdiff(required, names(table)), collapse = ", "))
  ranges <- GenomicRanges::GRanges(
    seqnames = table$chr,
    ranges = IRanges::IRanges(as.integer(table$start), as.integer(table$end))
  )
  ranges$gene <- as.character(table$gene)
  ranges
}

# Build mutually exclusive CRISPR experimental units. Both the legacy
 # FlowFISH tables and the EngreitzLab EPCrispr table can contain overlapping
 # intervals for the same gene. Connected, overlapping intervals are
 # collapsed within gene and chromosome. Mixed-label components cannot be
 # distinguished by a peak-level prediction and are excluded.
collapse_crispr_reference_units <- function(reference, drop_ambiguous = TRUE) {
  required <- c("chr", "start", "end", "gene", "label")
  if (!all(required %in% names(reference))) {
    stop("CRISPR reference lacks: ",
         paste(setdiff(required, names(reference)), collapse = ", "),
         call. = FALSE)
  }
  reference <- unique(reference[, required, drop = FALSE])
  reference$chr <- as.character(reference$chr)
  reference$gene <- as.character(reference$gene)
  reference$start <- as.integer(reference$start)
  reference$end <- as.integer(reference$end)
  reference$label <- as.character(reference$label)
  if (anyNA(reference[, c("chr", "gene", "start", "end", "label")])) {
    stop("CRISPR reference contains missing coordinates, gene, or label",
         call. = FALSE)
  }
  reference$group <- paste(reference$gene, reference$chr, sep = "\r")
  groups <- split(seq_len(nrow(reference)), reference$group)

  # Merge intervals with a small linear scan instead of repeatedly invoking
  # GRanges::reduce; this is substantially faster for the large EPCrispr table.
  collapsed_rows <- lapply(groups, function(index) {
    current <- reference[index, c("chr", "start", "end", "gene", "label"),
                          drop = FALSE]
    current <- current[order(current$start, current$end), , drop = FALSE]
    out <- list()
    out_index <- 1L
    current_start <- current$start[[1L]]
    current_end <- current$end[[1L]]
    current_labels <- current$label[[1L]]
    if (nrow(current) >= 2L) {
      for (i in seq.int(2L, nrow(current))) {
        if (current$start[[i]] <= current_end + 1L) {
          current_end <- max(current_end, current$end[[i]])
          current_labels <- c(current_labels, current$label[[i]])
        } else {
          out[[out_index]] <- data.frame(
            chr = current$chr[[1L]], start = current_start, end = current_end,
            gene = current$gene[[1L]],
            label = if (length(unique(current_labels)) == 1L) {
              current_labels[[1L]]
            } else "ambiguous",
            stringsAsFactors = FALSE
          )
          out_index <- out_index + 1L
          current_start <- current$start[[i]]
          current_end <- current$end[[i]]
          current_labels <- current$label[[i]]
        }
      }
    }
    out[[out_index]] <- data.frame(
      chr = current$chr[[1L]], start = current_start, end = current_end,
      gene = current$gene[[1L]],
      label = if (length(unique(current_labels)) == 1L) {
        current_labels[[1L]]
      } else "ambiguous",
      stringsAsFactors = FALSE
    )
    do.call(rbind, out)
  })
  collapsed <- do.call(rbind, collapsed_rows)
  units <- GenomicRanges::GRanges(
    seqnames = collapsed$chr,
    ranges = IRanges::IRanges(collapsed$start, collapsed$end),
    gene = collapsed$gene,
    label = collapsed$label
  )
  units$reference_id <- paste0("CRISPR_", seq_along(units))
  label_counts <- table(factor(as.character(units$label),
                               levels = c("positive", "negative", "ambiguous")))
  attr(units, "label_counts") <- label_counts
  if (drop_ambiguous) {
    units <- units[units$label != "ambiguous"]
    attr(units, "label_counts") <- label_counts
  }
  units
}

read_crispr_reference_units <- function(positive_path, negative_path,
                                        drop_ambiguous = TRUE) {
  read_labeled <- function(path, label) {
    x <- utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
    required <- c("chr", "start", "end", "gene")
    if (!all(required %in% names(x))) {
      stop("validation file lacks: ",
           paste(setdiff(required, names(x)), collapse = ", "), call. = FALSE)
    }
    x <- unique(x[, required, drop = FALSE])
    x$label <- label
    x
  }
  collapse_crispr_reference_units(rbind(
    read_labeled(positive_path, "positive"),
    read_labeled(negative_path, "negative")
  ), drop_ambiguous = drop_ambiguous)
}

# Cell-type mapping for the original Nasser et al. FlowFISH benchmark.  The
# raw table is hg19 and contains separate perturbation experiments for each
# cell line/condition; it should not be collapsed into one pooled reference
# before matching the dataset cell type.
flowfish_celltype_map <- list(
  Bcell = c("GM12878", "BJAB", "BJAB_anti-IgM_anti-CD40_4hr"),
  B = c("GM12878", "BJAB", "BJAB_anti-IgM_anti-CD40_4hr"),
  T = c("Jurkat", "Jurkat_anti-CD3_PMA_4hr"),
  # Keep the historical Tnk label usable for older PBMC result folders.  The
  # underlying experiments are Jurkat T cells; no NK-specific FlowFISH set is
  # present in this benchmark.
  Tnk = c("Jurkat", "Jurkat_anti-CD3_PMA_4hr"),
  CD4T = c("Jurkat", "Jurkat_anti-CD3_PMA_4hr"),
  Myeloid = c("THP1", "THP1_LPS_4hr"),
  # SHARE-seq LCL is represented by the GM12878 lymphoblastoid experiment.
  LCL = "GM12878",
  # K562 is retained as its own cell-line reference rather than being mixed
  # with THP1 myeloid experiments.
  K562 = "K562"
)

# Combine already-labelled CRISPR units and collapse overlapping intervals so
# the Fisher table has mutually exclusive experimental units.  This is used
# for the dataset-level `All` row after selecting only matched cell types.
combine_crispr_reference_units <- function(reference_list,
                                           drop_ambiguous = TRUE) {
  reference_list <- reference_list[!vapply(reference_list, is.null, logical(1L))]
  if (!length(reference_list)) return(NULL)
  frames <- lapply(reference_list, function(units) {
    if (!length(units)) return(NULL)
    data.frame(
      chr = as.character(GenomicRanges::seqnames(units)),
      start = as.integer(GenomicRanges::start(units)),
      end = as.integer(GenomicRanges::end(units)),
      gene = as.character(units$gene),
      label = as.character(units$label),
      stringsAsFactors = FALSE
    )
  })
  frames <- Filter(Negate(is.null), frames)
  if (!length(frames)) return(NULL)
  collapse_crispr_reference_units(
    do.call(rbind, frames), drop_ambiguous = drop_ambiguous
  )
}

# Read the original hg19 FlowFISH table and return one GRanges reference per
# broad dataset cell type.  Regulated==TRUE is the positive class and
# Regulated==FALSE is the tested negative class; this preserves the same label
# definition used to create the legacy positive/negative TSVs while retaining
# the original CellType information.
read_flowfish_reference_units <- function(path,
                                          celltype_map = flowfish_celltype_map,
                                          coordinate_chain = NULL,
                                          drop_ambiguous = TRUE) {
  if (!file.exists(path)) stop("Missing FlowFISH reference: ", path,
                               call. = FALSE)
  x <- data.table::fread(path, data.table = FALSE, comment.char = "#",
                         showProgress = FALSE)
  required <- c("chrPerturbationTarget", "startPerturbationTarget",
                "endPerturbationTarget", "GeneSymbol", "CellType",
                "Regulated")
  if (!all(required %in% names(x))) {
    stop("FlowFISH reference lacks: ",
         paste(setdiff(required, names(x)), collapse = ", "), call. = FALSE)
  }
  regulated <- toupper(as.character(x$Regulated))
  keep <- regulated %in% c("TRUE", "FALSE")
  x <- x[keep, , drop = FALSE]
  regulated <- toupper(as.character(x$Regulated))
  raw <- data.frame(
    chr = as.character(x$chrPerturbationTarget),
    start = as.integer(x$startPerturbationTarget),
    end = as.integer(x$endPerturbationTarget),
    gene = as.character(x$GeneSymbol),
    celltype = as.character(x$CellType),
    label = ifelse(regulated == "TRUE", "positive", "negative"),
    stringsAsFactors = FALSE
  )
  if (anyNA(raw[, c("chr", "start", "end", "gene", "celltype", "label")])) {
    stop("FlowFISH reference contains missing coordinates, gene, or celltype",
         call. = FALSE)
  }

  if (!is.null(coordinate_chain)) {
    if (is.character(coordinate_chain)) {
      coordinate_chain <- rtracklayer::import.chain(coordinate_chain)
    }
    source_ranges <- GenomicRanges::GRanges(
      seqnames = raw$chr,
      ranges = IRanges::IRanges(raw$start, raw$end),
      gene = raw$gene,
      celltype = raw$celltype,
      label = raw$label
    )
    lifted <- rtracklayer::liftOver(source_ranges, coordinate_chain)
    keep_lifted <- lengths(lifted) == 1L
    if (!any(keep_lifted)) stop("No FlowFISH intervals lifted to target build",
                               call. = FALSE)
    lifted_ranges <- unlist(lifted[keep_lifted], use.names = FALSE)
    lifted_meta <- S4Vectors::mcols(source_ranges[keep_lifted])
    raw <- data.frame(
      chr = as.character(GenomicRanges::seqnames(lifted_ranges)),
      start = as.integer(GenomicRanges::start(lifted_ranges)),
      end = as.integer(GenomicRanges::end(lifted_ranges)),
      gene = as.character(lifted_meta$gene),
      celltype = as.character(lifted_meta$celltype),
      label = as.character(lifted_meta$label),
      stringsAsFactors = FALSE
    )
  }

  references <- lapply(celltype_map, function(source_celltypes) {
    subset <- raw[raw$celltype %in% source_celltypes,
                  c("chr", "start", "end", "gene", "label"),
                  drop = FALSE]
    if (!nrow(subset)) return(NULL)
    collapse_crispr_reference_units(subset,
                                    drop_ambiguous = drop_ambiguous)
  })
  names(references) <- names(celltype_map)
  references[!vapply(references, is.null, logical(1L))]
}

# Read the GRCh38 enhancer-gene benchmark distributed with the EngreitzLab
# CRISPR comparison resource. Regulated supplies the binary experimental label;
# no SNP mapping is performed in this Fisher-test mode.
read_epcrispr_reference_units <- function(path,
                                          celltype = "K562",
                                          drop_ambiguous = TRUE) {
  if (!file.exists(path)) stop("Missing EPCrispr reference: ", path,
                               call. = FALSE)
  x <- data.table::fread(path, data.table = FALSE, showProgress = FALSE)
  required <- c("chrom", "chromStart", "chromEnd", "measuredGeneSymbol",
                "ValidConnection", "CellType", "Regulated")
  if (!all(required %in% names(x))) {
    stop("EPCrispr reference lacks: ",
         paste(setdiff(required, names(x)), collapse = ", "), call. = FALSE)
  }
  valid <- toupper(as.character(x$ValidConnection)) == "TRUE"
  cell_ok <- if (is.null(celltype) || !nzchar(celltype)) {
    rep(TRUE, nrow(x))
  } else {
    as.character(x$CellType) == celltype
  }
  regulated <- toupper(as.character(x$Regulated))
  keep <- valid & cell_ok & regulated %in% c("TRUE", "FALSE")
  x <- x[keep, , drop = FALSE]
  if (!nrow(x)) stop("No usable rows remain in EPCrispr reference", call. = FALSE)
  reference <- data.frame(
    chr = as.character(x$chrom),
    # EPCrispr chromStart/chromEnd are BED coordinates.
    start = as.integer(x$chromStart) + 1L,
    end = as.integer(x$chromEnd),
    gene = as.character(x$measuredGeneSymbol),
    label = ifelse(toupper(as.character(x$Regulated)) == "TRUE",
                   "positive", "negative"),
    stringsAsFactors = FALSE
  )
  collapse_crispr_reference_units(reference, drop_ambiguous = drop_ambiguous)
}

reference_units_hit <- function(pairs, reference_units, coordinate_chain = NULL) {
  if (!nrow(pairs) || !length(reference_units)) return(integer())
  query <- links_to_ranges(pairs, coordinate_chain = coordinate_chain)
  common <- base::intersect(
    GenomeInfoDb::seqlevelsInUse(query),
    GenomeInfoDb::seqlevelsInUse(reference_units)
  )
  if (!length(common)) return(integer())
  query <- GenomeInfoDb::keepSeqlevels(query, common, pruning.mode = "coarse")
  subject <- reference_units
  subject$unit_index_internal <- seq_along(subject)
  subject <- GenomeInfoDb::keepSeqlevels(subject, common,
                                         pruning.mode = "coarse")
  hits <- GenomicRanges::findOverlaps(query, subject, ignore.strand = TRUE)
  if (!length(hits)) return(integer())
  keep <- query$gene[S4Vectors::queryHits(hits)] ==
    subject$gene[S4Vectors::subjectHits(hits)]
  unique(subject$unit_index_internal[S4Vectors::subjectHits(hits)[keep]])
}

evaluate_crispr_reference_units <- function(universe, selected, reference_units,
                                            method, cutoff,
                                            benchmark = "CRISPR-FlowFISH",
                                            assessable_units = NULL,
                                            coordinate_chain = NULL) {
  universe <- unique(universe[, c("gene", "peak"), drop = FALSE])
  selected <- unique(selected[, c("gene", "peak"), drop = FALSE])
  universe_key <- paste(universe$gene, universe$peak, sep = "\r")
  selected_key <- paste(selected$gene, selected$peak, sep = "\r")
  selected <- universe[universe_key %in% selected_key, , drop = FALSE]

  assessable <- if (is.null(assessable_units)) {
    reference_units_hit(universe, reference_units, coordinate_chain)
  } else {
    unique(as.integer(assessable_units))
  }
  predicted <- base::intersect(
    reference_units_hit(selected, reference_units, coordinate_chain), assessable
  )
  labels <- reference_units$label
  positive <- assessable[labels[assessable] == "positive"]
  negative <- assessable[labels[assessable] == "negative"]

  tp <- length(base::intersect(predicted, positive))
  fp <- length(base::intersect(predicted, negative))
  fn <- length(base::setdiff(positive, predicted))
  tn <- length(base::setdiff(negative, predicted))
  contingency <- matrix(c(tp, fp, fn, tn), nrow = 2L, byrow = TRUE,
                        dimnames = list(
                          prediction = c("selected", "not_selected"),
                          validation = c("positive", "negative")
                        ))
  margins_defined <- all(rowSums(contingency) > 0) &&
    all(colSums(contingency) > 0)
  fisher <- if (margins_defined) {
    stats::fisher.test(contingency, alternative = "greater")
  } else NULL

  # Haldane-Anscombe is reported as a stable descriptive effect size.  The
  # exact one-sided Fisher p-value remains the inferential enrichment test.
  corrected <- contingency + 0.5
  if (margins_defined) {
    log_or_ha <- log(corrected[1, 1]) + log(corrected[2, 2]) -
      log(corrected[1, 2]) - log(corrected[2, 1])
    se_log_or_ha <- sqrt(sum(1 / corrected))
    or_ha <- exp(log_or_ha)
    ci_ha <- exp(log_or_ha + c(-1, 1) * stats::qnorm(0.975) * se_log_or_ha)
  } else {
    or_ha <- NA_real_
    ci_ha <- c(NA_real_, NA_real_)
  }
  precision <- if (tp + fp > 0L) tp / (tp + fp) else NA_real_
  recall <- if (tp + fn > 0L) tp / (tp + fn) else NA_real_
  specificity <- if (tn + fp > 0L) tn / (tn + fp) else NA_real_
  f1_denominator <- 2 * tp + fp + fn
  f1 <- if (f1_denominator > 0L) 2 * tp / f1_denominator else NA_real_
  counts <- attr(reference_units, "label_counts")

  data.frame(
    benchmark = benchmark,
    method = method,
    cutoff = cutoff,
    TP = tp, FP = fp, FN = fn, TN = tn,
    odds_ratio = if (margins_defined) unname(fisher$estimate) else NA_real_,
    CI_low = if (margins_defined) fisher$conf.int[[1L]] else NA_real_,
    CI_high = if (margins_defined) fisher$conf.int[[2L]] else NA_real_,
    pval = if (margins_defined) fisher$p.value else NA_real_,
    odds_ratio_ha = or_ha,
    CI_low_ha = ci_ha[[1L]],
    CI_high_ha = ci_ha[[2L]],
    precision = precision,
    recall = recall,
    specificity = specificity,
    F1 = f1,
    n_universe = nrow(universe),
    n_selected = nrow(selected),
    n_reference_units = length(reference_units),
    n_assessable_units = length(assessable),
    n_selected_units = length(predicted),
    n_positive_units = length(positive),
    n_negative_units = length(negative),
    n_ambiguous_units_excluded = if (is.null(counts)) NA_integer_ else
      unname(counts[["ambiguous"]]),
    stringsAsFactors = FALSE
  )
}

count_crispr_overlap <- function(query_pairs, validation_ranges, gene_match = TRUE) {
  if (!nrow(query_pairs) || !length(validation_ranges)) return(0L)
  query <- links_to_ranges(query_pairs)
  common <- intersect(
    GenomeInfoDb::seqlevelsInUse(query),
    GenomeInfoDb::seqlevelsInUse(validation_ranges)
  )
  if (!length(common)) return(0L)
  query <- GenomeInfoDb::keepSeqlevels(query, common, pruning.mode = "coarse")
  subject <- GenomeInfoDb::keepSeqlevels(validation_ranges, common, pruning.mode = "coarse")
  hits <- GenomicRanges::findOverlaps(query, subject, ignore.strand = TRUE)
  if (!length(hits)) return(0L)
  query_index <- S4Vectors::queryHits(hits)
  subject_index <- S4Vectors::subjectHits(hits)
  if (gene_match) {
    keep <- query$gene[query_index] == subject$gene[subject_index]
    query_index <- query_index[keep]
  }
  length(unique(query$pair_key[query_index]))
}

evaluate_binary_link_validation <- function(universe, selected, positive_ranges,
                                            negative_ranges, method, cutoff,
                                            benchmark = "CRISPR",
                                            positive_keys = NULL,
                                            negative_keys = NULL) {
  universe <- unique(universe[, c("gene", "peak"), drop = FALSE])
  selected <- unique(selected[, c("gene", "peak"), drop = FALSE])
  selected_key <- paste(selected$gene, selected$peak, sep = "\r")
  universe_key <- paste(universe$gene, universe$peak, sep = "\r")
  selected <- universe[paste(universe$gene, universe$peak, sep = "\r") %in% selected_key,
                       , drop = FALSE]
  non_selected <- universe[!universe_key %in% selected_key, , drop = FALSE]

  # Match the original incoming/validation implementation: count overlaps
  # separately on selected and non-selected tested pairs, without removing
  # pairs that happen to occur in both positive and negative references.
  tp <- count_crispr_overlap(selected, positive_ranges)
  fp <- count_crispr_overlap(selected, negative_ranges)
  fn <- count_crispr_overlap(non_selected, positive_ranges)
  tn <- count_crispr_overlap(non_selected, negative_ranges)

  contingency <- matrix(c(tp, fp, fn, tn), nrow = 2, byrow = TRUE,
  dimnames = list(prediction = c("selected", "not_selected"), 
                  validation = c("positive", "negative")))

fisher <- stats::fisher.test(contingency, alternative = "greater")
odds_defined <- all(rowSums(contingency) > 0) && all(colSums(contingency) > 0)

data.frame(
  benchmark = benchmark,
  method = method,
  cutoff = cutoff,
  TP = tp,
  FP = fp,
  FN = fn,
  TN = tn,
  odds_ratio = if (odds_defined) unname(fisher$estimate) else NA_real_,
  CI_low = if (odds_defined) fisher$conf.int[[1]] else NA_real_,
  CI_high = if (odds_defined) fisher$conf.int[[2]] else NA_real_,
  pval = fisher$p.value,
  n_universe = nrow(universe),
  n_selected = nrow(selected),
  stringsAsFactors = FALSE
)
}
