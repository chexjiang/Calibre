# PCHiC enrichment helpers for the standardized real-data benchmark outputs.
#
# The validation data are from Javierre et al. (Cell, 2016), DATA S1. The
# scMultiMap manuscript analysis uses the other-end (oe) fragment, the baitName
# gene, a cell-type-specific CHiCAGO score > 5, and a one-sided Fisher test.

pchic_method_specs <- data.frame(
  method = c("EpiLink-Poisson", "EpiLink-NB", "Poisson-GLM", "NB-GLM",
             "scMultiMap", "scMultiMap-noBatch", "Signac",
             "Signac (uncorrected)"),
  method_id = c("epilink_poisson", "epilink_nb", "poisson_glm", "nb_glm",
                "scmultimap", "scmultimap_nobatch", "signac",
                "signac_uncorrected"),
  stringsAsFactors = FALSE
)

# These mappings compare each broad EpiLink group with the closest available
# primary-cell PCHiC columns. FoeT is excluded because it is foetal thymus, and
# the Javierre matrix has no NK-specific column.
default_pchic_celltype_map <- list(
  Bcell = c("nB", "tB"),
  B = c("nB", "tB"),
  Myeloid = "Mon",
  Tnk = c("nCD4", "tCD4", "aCD4", "naCD4", "nCD8", "tCD8"),
  T = c("nCD4", "tCD4", "aCD4", "naCD4", "nCD8", "tCD8"),
  # NEAT-seq contains a CD4 T-cell population without a CD8/NK split.
  CD4T = c("nCD4", "tCD4", "aCD4", "naCD4")
)

# BMMC broad groups.  Javierre et al. has no dendritic-cell column, so
# Dendritic is intentionally omitted rather than assigned a misleading proxy.
# EP/Ery are the closest available erythroid progenitor/erythroid columns.
bmmc_pchic_celltype_map <- list(
  Bcell = c("nB", "tB"),
  Myeloid = c("Mon", "Mac0", "Mac1", "Mac2"),
  Erythroid = c("EP", "Ery"),
  T = c("nCD4", "tCD4", "aCD4", "naCD4", "nCD8", "tCD8")
)

# PBMC broad groups use T rather than the legacy Tnk label.
pbmc_pchic_celltype_map <- list(
  Bcell = c("nB", "tB"),
  Myeloid = "Mon",
  T = c("nCD4", "tCD4", "aCD4", "naCD4", "nCD8", "tCD8")
)

# Dogma stimulated/control datasets use the same broad Bcell/Myeloid/T labels
# as the PBMC benchmark, so use the same PCHiC mapping.
dogma_pchic_celltype_map <- pbmc_pchic_celltype_map

# SHARE-seq LCL is a lymphoblastoid/B-cell-line population.  The PCHiC
# reference has no LCL column, so nB/tB is an explicitly labelled primary
# B-cell proxy rather than an exact cell-line match.
share_pchic_celltype_map <- list(
  LCL = c("nB", "tB")
)

# Brain datasets do not have a matched Javierre PCHiC cell type.  Keep the
# validation runnable, but make the proxy explicit: each brain cell type uses
# the union of all available non-foetal PCHiC columns.  This is a broad
# genomic-reference enrichment, not a brain-cell-type-specific comparison.
brain_pchic_reference_columns <- c(
  "Mon", "Mac0", "Mac1", "Mac2", "Neu", "MK", "EP", "Ery",
  "nCD4", "tCD4", "aCD4", "naCD4", "nCD8", "tCD8", "nB", "tB"
)
brain_pchic_celltype_map <- setNames(
  rep(list(brain_pchic_reference_columns), 3L),
  c("Excitatory neurons", "Astrocytes", "Oligodendroglia")
)

load_pchic_reference <- function(project_root,
                                 celltype_map = default_pchic_celltype_map,
                                 chicago_cutoff = 5,
                                 chicago_operator = ">",
                                 expand_multigene_baits = FALSE) {
  data_path <- file.path(project_root, "reference", "PCHiC",
                         "PCHiC_peak_matrix_cutoff5.tsv")
  chain_path <- file.path(project_root, "reference", "PCHiC",
                          "hg19ToHg38.over.chain")
  if (!file.exists(data_path)) stop("Missing PCHiC matrix: ", data_path, call. = FALSE)
  if (!file.exists(chain_path)) stop("Missing liftOver chain: ", chain_path, call. = FALSE)
  if (!chicago_operator %in% c(">", ">=")) {
    stop("chicago_operator must be '>' or '>='", call. = FALSE)
  }

  score_columns <- unique(unlist(celltype_map, use.names = FALSE))
  columns <- c("oeChr", "oeStart", "oeEnd", "baitName", score_columns)
  pchic <- data.table::fread(data_path, select = columns, data.table = FALSE,
                            showProgress = FALSE)
  missing_columns <- setdiff(columns, names(pchic))
  if (length(missing_columns)) {
    stop("PCHiC matrix lacks: ", paste(missing_columns, collapse = ", "),
         call. = FALSE)
  }

  gr_hg19 <- GenomicRanges::makeGRangesFromDataFrame(
    pchic, seqnames.field = "oeChr", start.field = "oeStart",
    end.field = "oeEnd", keep.extra.columns = TRUE
  )
  GenomeInfoDb::seqlevelsStyle(gr_hg19) <- "UCSC"
  chain <- rtracklayer::import.chain(chain_path)
  lifted <- rtracklayer::liftOver(gr_hg19, chain)
  gr_hg38 <- unlist(lifted, use.names = FALSE)
  GenomeInfoDb::genome(gr_hg38) <- "hg38"

  make_reference <- function(score_names) {
    scores <- as.matrix(as.data.frame(S4Vectors::mcols(gr_hg38)[, score_names,
                                                                        drop = FALSE]))
    storage.mode(scores) <- "double"
    passes <- if (chicago_operator == ">") {
      rowSums(scores > chicago_cutoff, na.rm = TRUE) > 0L
    } else {
      rowSums(scores >= chicago_cutoff, na.rm = TRUE) > 0L
    }
    reference <- gr_hg38[passes]
    bait_names <- as.character(reference$baitName)
    if (expand_multigene_baits) {
      genes <- strsplit(bait_names, ";", fixed = TRUE)
      reference <- reference[rep(seq_along(reference), lengths(genes))]
      reference$gene <- unlist(genes, use.names = FALSE)
    } else {
      # Exact baitName matching reproduces scMultiMap's published code.
      reference$gene <- bait_names
    }
    reference
  }

  references <- lapply(celltype_map, make_reference)
  attr(references, "n_hg19") <- length(gr_hg19)
  attr(references, "n_hg38") <- length(gr_hg38)
  attr(references, "chicago_cutoff") <- chicago_cutoff
  attr(references, "chicago_operator") <- chicago_operator
  attr(references, "expand_multigene_baits") <- expand_multigene_baits
  references
}

read_pair_results <- function(path) {
  result <- data.table::fread(path, data.table = FALSE, showProgress = FALSE)
  required <- c("gene", "peak")
  if (!all(required %in% names(result))) {
    stop(path, " lacks: ", paste(setdiff(required, names(result)), collapse = ", "),
         call. = FALSE)
  }
  result
}

pair_key <- function(x) paste(x$gene, x$peak, sep = "\r")

restrict_pairs_to_universe <- function(pairs, universe) {
  if (!nrow(pairs)) return(universe[FALSE, c("gene", "peak"), drop = FALSE])
  keep <- pair_key(pairs) %in% pair_key(universe)
  unique(pairs[keep, c("gene", "peak"), drop = FALSE])
}

load_epilink_filter <- function(project_root) {
  package_environment <- new.env(parent = globalenv())
  source_files <- list.files(file.path(project_root, "R"), pattern = "[.]R$",
                             full.names = TRUE)
  for (source_file in source_files) sys.source(source_file, envir = package_environment)
  package_environment$EpiLinkFilter
}

select_method_pairs <- function(result, method_id, cutoff, universe,
                                method_dir, epilink_filter,
                                epilink_object = NULL, correction = TRUE) {
  if (grepl("^epilink_", method_id)) {
    object_path <- file.path(method_dir, "epilink_object.rds")
    if (!file.exists(object_path)) stop("Missing EpiLink object: ", object_path, call. = FALSE)
    object <- if (is.null(epilink_object)) readRDS(object_path) else epilink_object
    selected <- epilink_filter(
      object, stat = "zval", fdrcutoff = cutoff, correction = correction,
      conservative = FALSE, verbose = FALSE
    )
  } else if (method_id %in% validation_signac_method_ids) {
    signac <- restrict_pairs_to_universe(result, universe)
    result_index <- match(pair_key(signac), pair_key(result))
    signac$padj_common <- stats::p.adjust(
      validation_signac_pvalues(result, method_id)[result_index], method = "BH"
    )
    selected <- signac[is.finite(signac$padj_common) & signac$padj_common < cutoff,
                       c("gene", "peak"), drop = FALSE]
  } else {
    if (!"padj" %in% names(result)) stop(method_id, " result lacks padj", call. = FALSE)
    selected <- result[is.finite(result$padj) & result$padj < cutoff,
                       c("gene", "peak"), drop = FALSE]
  }
  restrict_pairs_to_universe(selected, universe)
}

positive_pchic_keys <- function(universe, reference, coordinate_chain = NULL) {
  if (!nrow(universe) || !length(reference)) return(character())
  query <- links_to_ranges(universe, coordinate_chain = coordinate_chain)
  common <- intersect(GenomeInfoDb::seqlevelsInUse(query),
                      GenomeInfoDb::seqlevelsInUse(reference))
  if (!length(common)) return(character())
  query <- GenomeInfoDb::keepSeqlevels(query, common, pruning.mode = "coarse")
  reference <- GenomeInfoDb::keepSeqlevels(reference, common, pruning.mode = "coarse")
  hits <- GenomicRanges::findOverlaps(query, reference, ignore.strand = TRUE)
  if (!length(hits)) return(character())
  query_index <- S4Vectors::queryHits(hits)
  reference_index <- S4Vectors::subjectHits(hits)
  matched <- query$gene[query_index] == reference$gene[reference_index]
  unique(query$pair_key[query_index[!is.na(matched) & matched]])
}

evaluate_pchic_enrichment <- function(universe, selected, reference, method,
                                      cutoff, dataset, celltype,
                                      positive_keys = NULL) {
  # scMultiMap's background contains tested pairs whose gene occurs in PCHiC.
  eligible <- universe$gene %in% unique(as.character(reference$gene))
  analysis_universe <- unique(universe[eligible, c("gene", "peak"), drop = FALSE])
  universe_keys <- pair_key(analysis_universe)
  selected_keys <- pair_key(selected)
  is_selected <- universe_keys %in% selected_keys
  positives <- if (is.null(positive_keys)) {
    positive_pchic_keys(analysis_universe, reference)
  } else positive_keys
  is_positive <- universe_keys %in% positives

  tp <- sum(is_selected & is_positive)
  fp <- sum(is_selected & !is_positive)
  fn <- sum(!is_selected & is_positive)
  tn <- sum(!is_selected & !is_positive)
  contingency <- matrix(c(tp, fp, fn, tn), nrow = 2L, byrow = TRUE,
                        dimnames = list(prediction = c("selected", "not_selected"),
                                        PCHiC = c("positive", "negative")))
  odds_defined <- all(rowSums(contingency) > 0L) && all(colSums(contingency) > 0L)
  fisher <- if (odds_defined) {
    stats::fisher.test(contingency, alternative = "greater")
  } else NULL

  data.frame(
    benchmark = "PCHiC",
    dataset = dataset,
    celltype = celltype,
    method = method,
    cutoff = cutoff,
    TP = tp, FP = fp, FN = fn, TN = tn,
    odds_ratio = if (odds_defined) unname(fisher$estimate) else NA_real_,
    CI_low = if (odds_defined) fisher$conf.int[[1L]] else NA_real_,
    CI_high = if (odds_defined) fisher$conf.int[[2L]] else NA_real_,
    pval = if (odds_defined) fisher$p.value else NA_real_,
    n_universe = nrow(analysis_universe),
    n_selected = sum(is_selected),
    n_pchic_links = sum(is_positive),
    stringsAsFactors = FALSE
  )
}

compute_dataset_pchic <- function(project_root, dataset, references,
                                  cutoffs = c(0.05, 0.10, 0.15, 0.20),
                                  allow_missing = TRUE, correction = TRUE,
                                  coordinate_chain = NULL,
                                  method_ids = pchic_method_specs$method_id) {
  dataset_root <- file.path(project_root, "results", "real", dataset)
  celltypes <- intersect(names(references),
                         list.dirs(dataset_root, recursive = FALSE, full.names = FALSE))
  if (!length(celltypes)) stop("No mapped cell-type results under ", dataset_root, call. = FALSE)
  epilink_filter <- load_epilink_filter(project_root)
  rows <- list()
  row_index <- 1L

  for (celltype in celltypes) {
    celltype_root <- file.path(dataset_root, celltype)
    universe_path <- file.path(celltype_root, "epilink_poisson", "all_pairs.tsv.gz")
    if (!file.exists(universe_path)) stop("Missing common universe: ", universe_path, call. = FALSE)
    universe <- unique(read_pair_results(universe_path)[, c("gene", "peak"), drop = FALSE])
    eligible_universe <- universe[
      universe$gene %in% unique(as.character(references[[celltype]]$gene)),
      c("gene", "peak"), drop = FALSE
    ]
    positive_keys <- positive_pchic_keys(
      eligible_universe, references[[celltype]], coordinate_chain = coordinate_chain
    )

    specs <- pchic_method_specs[pchic_method_specs$method_id %in% method_ids,
                                , drop = FALSE]
    for (method_index in seq_len(nrow(specs))) {
      method <- specs$method[method_index]
      method_id <- specs$method_id[method_index]
      method_dir <- file.path(celltype_root,
                              validation_result_method_id(method_id))
      result_path <- file.path(method_dir, "all_pairs.tsv.gz")
      if (!file.exists(result_path)) {
        if (allow_missing) next
        stop("Missing method result: ", result_path, call. = FALSE)
      }
      result <- read_pair_results(result_path)
      epilink_object <- if (grepl("^epilink_", method_id)) {
        readRDS(file.path(method_dir, "epilink_object.rds"))
      } else NULL
      for (cutoff in cutoffs) {
        selected <- select_method_pairs(result, method_id, cutoff, universe,
                                        method_dir, epilink_filter, epilink_object,
                                        correction = correction)
        rows[[row_index]] <- evaluate_pchic_enrichment(
          universe, selected, references[[celltype]], method, cutoff,
          dataset, celltype, positive_keys
        )
        row_index <- row_index + 1L
      }
    }
  }
  do.call(rbind, rows)
}
