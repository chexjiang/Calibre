# ArchR adapter for the model-based simulation benchmark.
#
# The simulation uses synthetic gene/peak identifiers (gene1, gene2, ...), so
# the adapter creates a small, deterministic synthetic genome annotation.  The
# count matrices and truth labels are unchanged; only genomic coordinates are
# supplied so that ArchR can construct its cis candidate set.

archr_patch_missing_imports <- local({
  patched <- FALSE
  function() {
    if (patched) return(invisible(FALSE))
    ns <- asNamespace("ArchR")
    # ArchR 1.0.3 currently calls SimpleList() in several namespace
    # functions but does not import it in NAMESPACE.  This is harmless in
    # installations where the import is present, but causes
    # "could not find function 'SimpleList'" on this R/Bioc installation.
    # Patch only the in-memory function bodies for this R session; the
    # installed ArchR package is not modified.
    function_names <- ls(ns, all.names = TRUE)
    for (function_name in function_names) {
      function_value <- tryCatch(
        get(function_name, envir = ns), error = function(e) NULL
      )
      if (!is.function(function_value)) next
      source_text <- paste(deparse(function_value), collapse = "\n")
      has_unqualified_simple_list <- grepl(
        "(?<![:[:alnum:]_.])SimpleList\\s*\\(", source_text,
        perl = TRUE
      )
      if (!has_unqualified_simple_list) next
      patched_text <- gsub(
        "(?<![:[:alnum:]_.])SimpleList\\s*\\(",
        "S4Vectors::SimpleList(", source_text, perl = TRUE
      )
      patched_function <- tryCatch(
        eval(parse(text = patched_text), envir = ns),
        error = function(e) NULL
      )
      if (is.function(patched_function)) {
        tryCatch(
          assignInNamespace(function_name, patched_function, ns = "ArchR"),
          error = function(e) NULL
        )
      }
    }
    patched <<- TRUE
    invisible(TRUE)
  }
})

archr_requirements <- function() {
  required <- c(
    "ArchR", "Matrix", "SingleCellExperiment", "SummarizedExperiment",
    "GenomicRanges", "GenomeInfoDb", "IRanges", "S4Vectors", "Rsamtools",
    "BSgenome.Hsapiens.UCSC.hg38"
  )
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      "ArchR simulation requires missing package(s): ",
      paste(missing, collapse = ", "), call. = FALSE
    )
  }
  archr_patch_missing_imports()
  invisible(TRUE)
}

archr_peak_ranges <- function(simulation) {
  peaks <- rownames(simulation$atac)
  pieces <- strsplit(peaks, "-", fixed = TRUE)
  if (any(lengths(pieces) != 3L)) {
    stop("Simulation peak names must have the form chr-start-end", call. = FALSE)
  }
  gr <- GenomicRanges::GRanges(
    seqnames = vapply(pieces, `[[`, character(1), 1L),
    ranges = IRanges::IRanges(
      start = as.integer(vapply(pieces, `[[`, character(1), 2L)),
      end = as.integer(vapply(pieces, `[[`, character(1), 3L))
    ),
    strand = "*"
  )
  names(gr) <- peaks
  S4Vectors::mcols(gr)$name <- peaks
  gr
}

archr_synthetic_annotations <- function(simulation, peak_ranges) {
  gene_ids <- rownames(simulation$rna)
  if (length(gene_ids) != length(peak_ranges)) {
    stop("Simulation RNA and ATAC feature counts do not match", call. = FALSE)
  }
  # Put each synthetic TSS at the center of its paired peak.  This preserves
  # the simulated pair universe while allowing ArchR's normal ±250-kb overlap
  # code to generate its own cis candidates.
  tss_position <- as.integer(round((
    GenomicRanges::start(peak_ranges) + GenomicRanges::end(peak_ranges)
  ) / 2))
  gene_ranges <- GenomicRanges::GRanges(
    seqnames = as.character(GenomicRanges::seqnames(peak_ranges)),
    ranges = IRanges::IRanges(
      start = pmax(1L, tss_position - 50L),
      end = tss_position + 50L
    ),
    strand = "*"
  )
  names(gene_ranges) <- gene_ids
  S4Vectors::mcols(gene_ranges)$symbol <- gene_ids
  exons <- gene_ranges
  tss <- GenomicRanges::GRanges(
    seqnames = as.character(GenomicRanges::seqnames(peak_ranges)),
    ranges = IRanges::IRanges(tss_position, width = 1L),
    strand = "*"
  )
  names(tss) <- gene_ids
  S4Vectors::mcols(tss)$symbol <- gene_ids
  list(
    gene_annotation = ArchR::createGeneAnnotation(
      genes = gene_ranges, exons = exons, TSS = tss
    ),
    gene_ranges = gene_ranges,
    tss = tss
  )
}

archr_synthetic_genome <- function(peak_ranges) {
  chromosomes <- unique(as.character(GenomicRanges::seqnames(peak_ranges)))
  ends <- vapply(chromosomes, function(chr) {
    max(GenomicRanges::end(
      peak_ranges[as.character(GenomicRanges::seqnames(peak_ranges)) == chr]
    ))
  }, numeric(1))
  chrom_sizes <- GenomicRanges::GRanges(
    seqnames = chromosomes,
    ranges = IRanges::IRanges(
      start = rep.int(1L, length(chromosomes)),
      end = as.integer(pmax(ends + 500000L, 1000000L))
    )
  )
  GenomeInfoDb::seqlengths(chrom_sizes) <- setNames(
    as.integer(pmax(ends + 500000L, 1000000L)), chromosomes
  )
  # ArchR uses the genome name to validate chromosome conventions.  Supplying
  # explicit chromosome sizes and an empty blacklist avoids downloading a
  # blacklist and keeps the synthetic coordinates deterministic.
  ArchR::createGenomeAnnotation(
    genome = "hg38", chromSizes = chrom_sizes,
    blacklist = GenomicRanges::GRanges(), filter = FALSE, filterChr = NULL
  )
}

archr_write_fragments <- function(simulation, peak_ranges, directory) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  occupied <- which(simulation$atac != 0, arr.ind = TRUE)
  if (!nrow(occupied)) {
    stop("Simulation has no accessible ATAC entries", call. = FALSE)
  }
  fragment_table <- data.frame(
    chr = as.character(GenomicRanges::seqnames(peak_ranges))[occupied[, 1]],
    start = GenomicRanges::start(peak_ranges)[occupied[, 1]],
    end = GenomicRanges::end(peak_ranges)[occupied[, 1]],
    barcode = colnames(simulation$atac)[occupied[, 2]],
    count = 1L,
    stringsAsFactors = FALSE
  )
  fragment_table <- fragment_table[order(
    fragment_table$chr, fragment_table$start, fragment_table$end,
    fragment_table$barcode
  ), , drop = FALSE]
  text_file <- file.path(directory, "sim.fragments.tsv")
  utils::write.table(
    fragment_table, text_file, sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = FALSE
  )
  gz_file <- Rsamtools::bgzip(text_file, overwrite = TRUE)
  Rsamtools::indexTabix(gz_file, format = "bed")
  gz_file
}

archr_make_sce <- function(simulation, gene_ranges) {
  counts <- Matrix::Matrix(simulation$rna, sparse = TRUE)
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts), rowRanges = gene_ranges
  )
  colnames(sce) <- paste0("sim#", colnames(simulation$rna))
  rownames(sce) <- rownames(simulation$rna)
  SummarizedExperiment::rowData(sce)$name <- rownames(simulation$rna)
  sce
}

archr_extract_selected <- function(project, simulation, q, cor_cutoff = 0.45) {
  links <- ArchR::getPeak2GeneLinks(
    project, corCutOff = cor_cutoff, FDRCutOff = q,
    varCutOffATAC = NULL, varCutOffRNA = NULL,
    returnLoops = FALSE
  )
  if (is.null(links) || !nrow(links)) {
    return(list(selected = integer(), n_fitted = 0L))
  }
  peak_ranges <- ArchR::getPeakSet(project)
  gene_ranges <- ArchR::getGenes(project)
  peak_names <- names(peak_ranges)
  if ("name" %in% names(S4Vectors::mcols(peak_ranges))) {
    stored_names <- as.character(S4Vectors::mcols(peak_ranges)$name)
    if (all(stored_names %in% rownames(simulation$atac))) {
      peak_names <- stored_names
    }
  }
  if (is.null(peak_names)) peak_names <- rownames(simulation$atac)
  gene_symbols <- as.character(S4Vectors::mcols(gene_ranges)$symbol)
  if (is.null(gene_symbols) || !length(gene_symbols)) {
    gene_symbols <- rownames(simulation$rna)
  }
  linked_peak <- peak_names[as.integer(links$idxATAC)]
  linked_gene <- gene_symbols[as.integer(links$idxRNA)]
  candidate_key <- paste(
    rownames(simulation$rna), rownames(simulation$atac), sep = "\r"
  )
  linked_key <- paste(linked_gene, linked_peak, sep = "\r")
  list(
    selected = which(candidate_key %in% linked_key),
    # The ArchR object tests all cis pairs within maxDist.  The benchmark
    # score is restricted to the fixed candidate universe, as for other methods.
    n_fitted = length(unique(linked_key))
  )
}

#' Run ArchR Peak2GeneLinks on one model-based simulation dataset.
#'
#' This intentionally uses the direct GeneExpressionMatrix because the
#' simulation supplies matched RNA and ATAC cells.  It does not alter the
#' simulated counts, truth labels, or candidate pair universe.
run_archr_simulation <- function(simulation, q = 0.1, seed = 1L,
                                  n_cores = 1L, work_directory = NULL,
                                  cor_cutoff = 0.45, max_dist = 250000L,
                                  k = 100L) {
  archr_requirements()
  if (!is.null(work_directory)) {
    dir.create(work_directory, recursive = TRUE, showWarnings = FALSE)
  } else {
    work_directory <- tempfile("epilink-archr-")
    dir.create(work_directory, recursive = TRUE, showWarnings = FALSE)
  }
  on.exit(unlink(work_directory, recursive = TRUE, force = TRUE), add = TRUE)
  started <- proc.time()[["elapsed"]]
  peak_ranges <- archr_peak_ranges(simulation)
  annotations <- archr_synthetic_annotations(simulation, peak_ranges)
  genome_annotation <- archr_synthetic_genome(peak_ranges)
  fragments <- archr_write_fragments(
    simulation, peak_ranges, file.path(work_directory, "fragments")
  )
  arrow_directory <- file.path(work_directory, "arrow")
  project_directory <- file.path(work_directory, "project")
  dir.create(arrow_directory, recursive = TRUE, showWarnings = FALSE)
  dir.create(project_directory, recursive = TRUE, showWarnings = FALSE)

  ArchR::addArchRThreads(threads = as.integer(n_cores))
  ArchR::addArchRLocking(locking = TRUE)
  arrow_files <- ArchR::createArrowFiles(
    inputFiles = fragments, sampleNames = "sim",
    outputNames = file.path(arrow_directory, "sim"),
    geneAnnotation = annotations$gene_annotation,
    genomeAnnotation = genome_annotation,
    minTSS = 0, minFrags = 0, maxFrags = 1e9,
    addTileMat = FALSE, addGeneScoreMat = FALSE,
    filterFrags = 0, filterTSS = 0,
    threads = as.integer(n_cores), subThreading = FALSE,
    force = TRUE, verbose = FALSE,
    logFile = file.path(work_directory, "createArrowFiles.log")
  )
  project <- ArchR::ArchRProject(
    ArrowFiles = arrow_files,
    outputDirectory = project_directory,
    copyArrows = FALSE,
    geneAnnotation = annotations$gene_annotation,
    genomeAnnotation = genome_annotation,
    showLogo = FALSE
  )
  project <- ArchR::addPeakSet(project, peakSet = peak_ranges, force = TRUE)
  project <- ArchR::addPeakMatrix(
    project, force = TRUE, threads = as.integer(n_cores),
    verbose = FALSE,
    logFile = file.path(work_directory, "addPeakMatrix.log")
  )
  project <- ArchR::addGeneExpressionMatrix(
    input = project, seRNA = archr_make_sce(simulation, annotations$gene_ranges),
    strictMatch = TRUE, scaleTo = 10000, force = TRUE,
    threads = as.integer(n_cores), verbose = FALSE,
    logFile = file.path(work_directory, "addGeneExpressionMatrix.log")
  )
  dims <- max(2L, min(10L, simulation$n_cells - 1L, simulation$n_pairs - 1L))
  project <- ArchR::addIterativeLSI(
    ArchRProj = project, useMatrix = "PeakMatrix",
    name = "IterativeLSI", iterations = 1L,
    varFeatures = max(2L, min(250L, simulation$n_pairs)),
    dimsToUse = seq_len(dims),
    outlierQuantiles = NULL, filterBias = FALSE,
    sampleCellsPre = NULL, sampleCellsFinal = NULL, keep0lsi = TRUE,
    clusterParams = list(
      resolution = 0.2, sampleCells = simulation$n_cells, n.start = 1L
    ),
    threads = as.integer(n_cores), seed = as.integer(seed),
    force = TRUE, verbose = FALSE,
    logFile = file.path(work_directory, "addIterativeLSI.log")
  )
  # ArchR's paper uses k=100 for very large datasets.  The present benchmark
  # includes only 100--400 cells, so use at most one quarter of the cells to
  # retain multiple KNN groups (while never exceeding the paper's k=100).
  group_k <- as.integer(max(2L, min(k, floor(simulation$n_cells / 4L),
                                     simulation$n_cells - 1L)))
  project <- ArchR::addPeak2GeneLinks(
    ArchRProj = project, reducedDims = "IterativeLSI",
    useMatrix = "GeneExpressionMatrix", dimsToUse = seq_len(dims),
    k = group_k, knnIteration = min(500L, simulation$n_cells),
    overlapCutoff = 0.8, maxDist = as.integer(max_dist),
    scaleTo = 10000, log2Norm = TRUE, predictionCutoff = 0.0,
    addEmpiricalPval = FALSE, addPermutedPval = FALSE,
    seed = as.integer(seed), threads = as.integer(n_cores),
    verbose = FALSE,
    logFile = file.path(work_directory, "addPeak2GeneLinks.log")
  )
  selected <- archr_extract_selected(project, simulation, q, cor_cutoff = cor_cutoff)
  score <- score_simulation_selection(selected$selected, simulation$truth)
  score$method <- "ArchR"
  score$method_id <- "archr"
  score$finite_sample_correction <- FALSE
  score$plus_one <- FALSE
  score$n_fitted <- selected$n_fitted
  score$runtime_seconds <- proc.time()[["elapsed"]] - started
  score$error <- NA_character_
  attr(score, "selected") <- selected$selected
  score
}
