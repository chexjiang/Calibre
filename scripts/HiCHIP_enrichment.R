# H3K27ac HiChIP enrichment helpers following the SCENT consistency analysis.
# Coordinates in the workbook are hg19 BED-like intervals. They are lifted to
# hg38, promoters are +/- 1 kb around gene TSSs, and only enhancer--promoter
# loops (one promoter anchor and one non-promoter anchor) are retained.

hichip_method_specs <- data.frame(
  method = c("EpiLink-Poisson", "EpiLink-NB", "Poisson-GLM", "NB-GLM",
             "scMultiMap", "scMultiMap-noBatch", "Signac",
             "Signac (uncorrected)"),
  method_id = c("epilink_poisson", "epilink_nb", "poisson_glm", "nb_glm",
                "scmultimap", "scmultimap_nobatch", "signac",
                "signac_uncorrected"), stringsAsFactors = FALSE
)

default_hichip_sheet_map <- list(
  Tnk = c("Naive H3K27ac Loops", "Th17 H3K27ac Loops", "Treg H3K27ac Loops"),
  T = c("Naive H3K27ac Loops", "Th17 H3K27ac Loops", "Treg H3K27ac Loops"),
  Bcell = "GM12878 H3K27ac Loops",
  B = "GM12878 H3K27ac Loops",
  # NEAT-seq has a broad CD4T population; pool the three CD4 T-cell
  CD4T = c("Naive H3K27ac Loops", "Th17 H3K27ac Loops", "Treg H3K27ac Loops"),
  K562 = "K562 H3K27ac Loops"
)

# BMMC broad groups supported by the available workbook.  There are no
# matched Myeloid, Dendritic, or Erythroid sheets in this HiChIP reference.
bmmc_hichip_sheet_map <- list(
  Bcell = "GM12878 H3K27ac Loops",
  T = c("Naive H3K27ac Loops", "Th17 H3K27ac Loops", "Treg H3K27ac Loops")
)

# PBMC broad groups use T rather than the legacy Tnk label.
pbmc_hichip_sheet_map <- list(
  Bcell = "GM12878 H3K27ac Loops",
  T = c("Naive H3K27ac Loops", "Th17 H3K27ac Loops", "Treg H3K27ac Loops")
)

# Dogma stimulated/control datasets use the same broad Bcell/Myeloid/T labels
# as PBMC.  No matched Myeloid sheet exists, so Myeloid uses the released-sheet
# union as an explicit broad proxy rather than being silently dropped.
dogma_hichip_sheet_map <- c(
  pbmc_hichip_sheet_map,
  list(Myeloid = c(
    "mES H3K27ac Loops", "Naive H3K27ac Loops", "Th17 H3K27ac Loops",
    "Treg H3K27ac Loops", "GM12878 H3K27ac Loops", "K562 H3K27ac Loops",
    "My-La H3K27ac Loops", "HCASMC H3K27ac Loops"
  ))
)

# SHARE-seq LCL is matched to the GM12878 H3K27ac HiChIP sheet.
share_hichip_sheet_map <- list(
  LCL = "GM12878 H3K27ac Loops"
)

# No brain-matched H3K27ac HiChIP sheet is present in the SCENT workbook.
# Use all released sheets as an explicitly broad proxy for each Brain cell
# type, while retaining the original Brain labels in the output table.
brain_hichip_sheets <- c(
  "mES H3K27ac Loops", "Naive H3K27ac Loops", "Th17 H3K27ac Loops",
  "Treg H3K27ac Loops", "GM12878 H3K27ac Loops", "K562 H3K27ac Loops",
  "My-La H3K27ac Loops", "HCASMC H3K27ac Loops"
)
brain_hichip_sheet_map <- setNames(
  rep(list(brain_hichip_sheets), 3L),
  c("Excitatory neurons", "Astrocytes", "Oligodendroglia")
)

hichip_pair_key <- function(x) paste(x$gene, x$peak, sep = "\r")

hichip_restrict_pairs <- function(pairs, universe) {
  if (!nrow(pairs)) return(universe[FALSE, c("gene", "peak"), drop = FALSE])
  keep <- hichip_pair_key(pairs) %in% hichip_pair_key(universe)
  unique(pairs[keep, c("gene", "peak"), drop = FALSE])
}

hichip_read_pairs <- function(path) {
  x <- data.table::fread(path, data.table = FALSE, showProgress = FALSE)
  if (!all(c("gene", "peak") %in% names(x))) {
    stop(path, " lacks gene and/or peak columns", call. = FALSE)
  }
  x
}

hichip_make_promoters <- function() {
  if (!requireNamespace("EnsDb.Hsapiens.v86", quietly = TRUE) ||
      !requireNamespace("ensembldb", quietly = TRUE)) {
    stop("EnsDb.Hsapiens.v86 and ensembldb are required for HiChIP TSS annotation",
         call. = FALSE)
  }
  genes <- ensembldb::genes(
    EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86,
    return.type = "GRanges",
    columns = c("gene_id", "gene_name", "gene_biotype")
  )
  levels <- GenomeInfoDb::seqlevels(genes)
  GenomeInfoDb::seqlevels(genes) <- ifelse(levels == "MT", "chrM", paste0("chr", levels))
  if (!"gene_name" %in% names(S4Vectors::mcols(genes))) {
    stop("EnsDb annotation lacks gene_name", call. = FALSE)
  }
  keep <- !is.na(genes$gene_name) & nzchar(as.character(genes$gene_name))
  genes <- genes[keep & !duplicated(as.character(genes$gene_id))]
  promoters <- suppressWarnings(GenomicRanges::promoters(
    genes, upstream = 1000L, downstream = 1001L
  ))
  promoters$gene <- as.character(genes$gene_name)
  promoters
}

hichip_read_sheet <- function(path, sheet) {
  x <- as.data.frame(suppressWarnings(readxl::read_excel(path, sheet = sheet)),
                     stringsAsFactors = FALSE)
  required <- c("chr_x", "start_x", "stop_x", "chr_y", "start_y", "stop_y")
  if (!all(required %in% names(x))) stop("HiChIP sheet lacks required coordinates: ", sheet)
  x <- x[stats::complete.cases(x[, required, drop = FALSE]), required, drop = FALSE]
  x$chr_x <- ifelse(grepl("^chr", x$chr_x), as.character(x$chr_x), paste0("chr", x$chr_x))
  x$chr_y <- ifelse(grepl("^chr", x$chr_y), as.character(x$chr_y), paste0("chr", x$chr_y))
  # Workbook intervals are BED-like (0-based start, half-open end).
  x$start_x <- as.integer(x$start_x) + 1L; x$stop_x <- as.integer(x$stop_x)
  x$start_y <- as.integer(x$start_y) + 1L; x$stop_y <- as.integer(x$stop_y)
  x <- x[x$start_x <= x$stop_x & x$start_y <= x$stop_y, , drop = FALSE]
  x$loop_id <- seq_len(nrow(x)); x
}

hichip_lift_loops <- function(loops, chain) {
  anchors <- rbind(
    data.frame(loop_id = loops$loop_id, side = "x", chr = loops$chr_x,
               start = loops$start_x, end = loops$stop_x),
    data.frame(loop_id = loops$loop_id, side = "y", chr = loops$chr_y,
               start = loops$start_y, end = loops$stop_y)
  )
  gr <- GenomicRanges::GRanges(anchors$chr, IRanges::IRanges(anchors$start, anchors$end))
  GenomeInfoDb::seqlevelsStyle(gr) <- "UCSC"
  mapped <- rtracklayer::liftOver(gr, chain)
  keep <- lengths(mapped) == 1L
  if (!any(keep)) return(anchors[FALSE, , drop = FALSE])
  lifted <- unlist(mapped[keep], use.names = FALSE)
  out <- anchors[keep, , drop = FALSE]
  out$chr <- as.character(GenomicRanges::seqnames(lifted))
  out$start <- GenomicRanges::start(lifted); out$end <- GenomicRanges::end(lifted)
  out
}

hichip_make_reference <- function(path, promoters,
                                  sheet_map = default_hichip_sheet_map) {
  chain_path <- file.path(dirname(path), "hg19ToHg38.over.chain")
  if (!file.exists(path) || !file.exists(chain_path)) {
    stop("HiChIP workbook or hg19ToHg38.over.chain is missing", call. = FALSE)
  }
  chain <- rtracklayer::import.chain(chain_path)
  refs <- list()
  for (celltype in names(sheet_map)) {
    link_parts <- list()
    for (sheet in sheet_map[[celltype]]) {
      loops <- hichip_read_sheet(path, sheet)
      anchors <- hichip_lift_loops(loops, chain)
      if (!nrow(anchors)) next
      anchor_gr <- GenomicRanges::GRanges(
        anchors$chr, IRanges::IRanges(anchors$start, anchors$end)
      )
      hits <- GenomicRanges::findOverlaps(anchor_gr, promoters, ignore.strand = TRUE)
      genes_by_key <- split(
        as.character(promoters$gene[S4Vectors::subjectHits(hits)]),
        paste(anchors$loop_id[S4Vectors::queryHits(hits)],
              anchors$side[S4Vectors::queryHits(hits)], sep = ":")
      )
      for (indices in split(seq_len(nrow(anchors)), anchors$loop_id)) {
        xi <- indices[anchors$side[indices] == "x"]
        yi <- indices[anchors$side[indices] == "y"]
        if (length(xi) != 1L || length(yi) != 1L) next
        xkey <- paste(anchors$loop_id[xi], "x", sep = ":")
        ykey <- paste(anchors$loop_id[yi], "y", sep = ":")
        xgenes <- unique(as.character(genes_by_key[[xkey]]))
        ygenes <- unique(as.character(genes_by_key[[ykey]]))
        xgenes <- xgenes[!is.na(xgenes) & nzchar(xgenes)]
        ygenes <- ygenes[!is.na(ygenes) & nzchar(ygenes)]
        if (length(xgenes) && !length(ygenes)) {
          for (gene in xgenes) link_parts[[length(link_parts) + 1L]] <-
            data.frame(gene = gene, chr = anchors$chr[yi], start = anchors$start[yi],
                       end = anchors$end[yi], stringsAsFactors = FALSE)
        } else if (length(ygenes) && !length(xgenes)) {
          for (gene in ygenes) link_parts[[length(link_parts) + 1L]] <-
            data.frame(gene = gene, chr = anchors$chr[xi], start = anchors$start[xi],
                       end = anchors$end[xi], stringsAsFactors = FALSE)
        }
      }
    }
    if (!length(link_parts)) {
      refs[[celltype]] <- GenomicRanges::GRanges()
    } else {
      links <- unique(do.call(rbind, link_parts))
      gr <- GenomicRanges::GRanges(links$chr, IRanges::IRanges(links$start, links$end))
      gr$gene <- links$gene; refs[[celltype]] <- gr
    }
  }
  refs
}

hichip_positive_keys <- function(universe, reference, coordinate_chain = NULL) {
  if (!nrow(universe) || !length(reference)) return(character())
  query <- links_to_ranges(universe, coordinate_chain = coordinate_chain)
  common <- intersect(GenomeInfoDb::seqlevelsInUse(query), GenomeInfoDb::seqlevelsInUse(reference))
  if (!length(common)) return(character())
  query <- GenomeInfoDb::keepSeqlevels(query, common, pruning.mode = "coarse")
  reference <- GenomeInfoDb::keepSeqlevels(reference, common, pruning.mode = "coarse")
  hits <- GenomicRanges::findOverlaps(query, reference, ignore.strand = TRUE)
  if (!length(hits)) return(character())
  qi <- S4Vectors::queryHits(hits); si <- S4Vectors::subjectHits(hits)
  unique(query$pair_key[query$gene[qi] == reference$gene[si]])
}

hichip_select_pairs <- function(result, method_id, cutoff, universe, method_dir,
                                epilink_filter, correction, epilink_object = NULL) {
  if (grepl("^epilink_", method_id)) {
    object <- if (is.null(epilink_object)) readRDS(file.path(method_dir, "epilink_object.rds")) else epilink_object
    selected <- epilink_filter(object, stat = "zval", fdrcutoff = cutoff,
                               correction = correction, conservative = FALSE, verbose = FALSE)
  } else if (method_id %in% validation_signac_method_ids) {
    selected <- hichip_restrict_pairs(result, universe)
    if (!nrow(selected)) return(selected)
    idx <- match(hichip_pair_key(selected), hichip_pair_key(result))
    selected$padj_hichip <- stats::p.adjust(
      validation_signac_pvalues(result, method_id)[idx], method = "BH"
    )
    selected <- selected[is.finite(selected$padj_hichip) & selected$padj_hichip < cutoff,
                         c("gene", "peak"), drop = FALSE]
  } else {
    selected <- result[is.finite(result$padj) & result$padj < cutoff,
                       c("gene", "peak"), drop = FALSE]
  }
  hichip_restrict_pairs(selected, universe)
}

hichip_evaluate <- function(universe, selected, reference, method, cutoff,
                            dataset, celltype, positive_keys) {
  eligible <- universe$gene %in% unique(as.character(reference$gene))
  test_universe <- unique(universe[eligible, c("gene", "peak"), drop = FALSE])
  keys <- hichip_pair_key(test_universe)
  selected_keys <- hichip_pair_key(selected)
  is_selected <- keys %in% selected_keys; is_positive <- keys %in% positive_keys
  tp <- sum(is_selected & is_positive); fp <- sum(is_selected & !is_positive)
  fn <- sum(!is_selected & is_positive); tn <- sum(!is_selected & !is_positive)
  tab <- matrix(c(tp, fp, fn, tn), nrow = 2L, byrow = TRUE)
  defined <- all(rowSums(tab) > 0L) && all(colSums(tab) > 0L)
  test <- if (defined) stats::fisher.test(tab, alternative = "greater") else NULL
  data.frame(benchmark = "HiChIP", dataset = dataset, celltype = celltype,
             method = method, cutoff = cutoff, TP = tp, FP = fp, FN = fn, TN = tn,
             odds_ratio = if (defined) unname(test$estimate) else NA_real_,
             CI_low = if (defined) test$conf.int[[1L] ] else NA_real_,
             CI_high = if (defined) test$conf.int[[2L] ] else NA_real_,
             pval = if (defined) test$p.value else NA_real_,
             n_universe = nrow(test_universe), n_selected = sum(is_selected),
             n_hichip_links = sum(is_positive), stringsAsFactors = FALSE)
}

compute_dataset_hichip <- function(project_root, dataset, references,
                                   cutoffs = c(0.05, 0.10, 0.15, 0.20),
                                   correction = TRUE, allow_missing = TRUE,
                                   coordinate_chain = NULL,
                                   method_ids = hichip_method_specs$method_id) {
  root <- file.path(project_root, "results", "real", dataset)
  celltypes <- intersect(names(references), list.dirs(root, recursive = FALSE, full.names = FALSE))
  if (!length(celltypes)) stop("No mapped cell-type results under ", root, call. = FALSE)
  env <- new.env(parent = globalenv())
  for (f in list.files(file.path(project_root, "R"), pattern = "[.]R$", full.names = TRUE)) sys.source(f, envir = env)
  epilink_filter <- env$EpiLinkFilter
  rows <- list(); row_i <- 1L
  for (celltype in celltypes) {
    cell_root <- file.path(root, celltype)
    universe_path <- file.path(cell_root, "epilink_poisson", "all_pairs.tsv.gz")
    if (!file.exists(universe_path)) stop("Missing common universe: ", universe_path, call. = FALSE)
    universe <- unique(hichip_read_pairs(universe_path)[, c("gene", "peak"), drop = FALSE])
    positives <- hichip_positive_keys(
      universe, references[[celltype]], coordinate_chain = coordinate_chain
    )
    specs <- hichip_method_specs[hichip_method_specs$method_id %in% method_ids,
                                 , drop = FALSE]
    for (i in seq_len(nrow(specs))) {
      method <- specs$method[i]; method_id <- specs$method_id[i]
      method_dir <- file.path(cell_root, validation_result_method_id(method_id)); result_path <- file.path(method_dir, "all_pairs.tsv.gz")
      if (!file.exists(result_path)) {
        if (allow_missing) next
        stop("Missing method result: ", result_path, call. = FALSE)
      }
      result <- hichip_read_pairs(result_path)
      object <- if (grepl("^epilink_", method_id)) readRDS(file.path(method_dir, "epilink_object.rds")) else NULL
      for (cutoff in cutoffs) {
        selected <- hichip_select_pairs(result, method_id, cutoff, universe, method_dir, epilink_filter, correction, object)
        rows[[row_i]] <- hichip_evaluate(universe, selected, references[[celltype]], method, cutoff, dataset, celltype, positives)
        row_i <- row_i + 1L
      }
    }
  }
  if (!length(rows)) stop("No HiChIP enrichment rows were generated", call. = FALSE)
  do.call(rbind, rows)
}
