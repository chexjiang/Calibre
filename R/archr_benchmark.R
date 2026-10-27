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
    # Some ArchR releases call symbols that are not imported into the
    # namespace.  This is harmless in installations where the imports are
    # present, but fails at runtime on otherwise valid Bioconductor setups.
    # Keep the compatibility patch in-memory; never modify the installed
    # ArchR package.
    replacement_map <- c(
      SimpleList = "S4Vectors::SimpleList",
      h5testFileLocking = "rhdf5::h5testFileLocking",
      h5disableFileLocking = "rhdf5::h5disableFileLocking",
      h5enableFileLocking = "rhdf5::h5enableFileLocking",
      h5closeAll = "rhdf5::h5closeAll",
      h5read = "rhdf5::h5read",
      h5write = "rhdf5::h5write",
      h5createFile = "rhdf5::h5createFile",
      h5createDataset = "rhdf5::h5createDataset",
      h5writeDataset = "rhdf5::h5writeDataset",
      h5createGroup = "rhdf5::h5createGroup",
      h5delete = "rhdf5::h5delete",
      h5ls = "rhdf5::h5ls",
      h5readAttributes = "rhdf5::h5readAttributes",
      h5writeAttribute = "rhdf5::h5writeAttribute",
      H5Fopen = "rhdf5::H5Fopen",
      H5Fclose = "rhdf5::H5Fclose",
      H5Fcreate = "rhdf5::H5Fcreate",
      H5Lcreate_external = "rhdf5::H5Lcreate_external",
      H5Pcreate = "rhdf5::H5Pcreate",
      H5Dopen = "rhdf5::H5Dopen",
      H5Dread = "rhdf5::H5Dread",
      H5Dclose = "rhdf5::H5Dclose",
      SummarizedExperiment = "SummarizedExperiment::SummarizedExperiment",
      assay = "SummarizedExperiment::assay",
      rowRanges = "SummarizedExperiment::rowRanges",
      rowData = "SummarizedExperiment::rowData",
      colData = "SummarizedExperiment::colData",
      assays = "SummarizedExperiment::assays"
    )
    # Patch only the in-memory function bodies for this R session; the
    # installed ArchR package is not modified.
    function_names <- ls(ns, all.names = TRUE)
    for (function_name in function_names) {
      function_value <- tryCatch(
        get(function_name, envir = ns), error = function(e) NULL
      )
      if (!is.function(function_value)) next
      source_text <- paste(deparse(function_value), collapse = "\n")
      has_unqualified_pipe <- grepl(
        "(?<![:[:alnum:]_.`])%>%", source_text, perl = TRUE
      )
      patched_text <- source_text
      for (symbol in names(replacement_map)) {
        patched_text <- gsub(
          paste0("(?<![:[:alnum:]_.])", symbol, "\\s*\\("),
          paste0(replacement_map[[symbol]], "("),
          patched_text, perl = TRUE
        )
      }
      if (identical(patched_text, source_text) && !has_unqualified_pipe) next
      patched_function <- tryCatch(
        eval(parse(text = patched_text), envir = ns),
        error = function(e) NULL
      )
      if (is.function(patched_function)) {
        if (has_unqualified_pipe) {
          # Namespace environments are locked, so provide the missing
          # operator through the function's lexical environment instead.
          patched_environment <- new.env(parent = environment(patched_function))
          assign(
            "%>%", get("%>%", envir = asNamespace("magrittr")),
            envir = patched_environment
          )
          environment(patched_function) <- patched_environment
        }
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
    "rhdf5", "magrittr", "BSgenome.Hsapiens.UCSC.hg38"
  )
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      "ArchR runs require missing package(s): ",
      paste(missing, collapse = ", "), call. = FALSE
    )
  }
  archr_patch_missing_imports()
  invisible(TRUE)
}

archr_var_features <- function(n_features, cap = 5000L) {
  n_features <- as.integer(n_features)
  if (is.na(n_features) || n_features <= 1000L) {
    stop(
      "ArchR iterative LSI requires more than 1000 peak features; found ",
      n_features, call. = FALSE
    )
  }
  min(as.integer(cap), n_features)
}

archr_lsi_dims <- function(n_cells, n_features, cap = 30L) {
  n_cells <- as.integer(n_cells)
  n_features <- as.integer(n_features)
  if (is.na(n_cells) || n_cells <= 4L) {
    stop("ArchR iterative LSI requires more than four cells", call. = FALSE)
  }
  max(2L, min(as.integer(cap), n_cells - 3L, n_features - 1L))
}

archr_cell_matches_sample <- function(cell, sample) {
  cell <- as.character(cell)
  sample <- as.character(sample)
  startsWith(cell, paste0(sample, "_")) ||
    startsWith(cell, paste0(sample, "#")) ||
    endsWith(cell, paste0("-", sample))
}

archr_strip_sample_from_cell <- function(cell, sample) {
  cell <- as.character(cell)
  sample <- as.character(sample)
  prefix <- paste0(sample, "_")
  hash_prefix <- paste0(sample, "#")
  suffix <- paste0("-", sample)
  if (startsWith(cell, prefix)) {
    return(sub(prefix, "", cell, fixed = TRUE))
  }
  if (startsWith(cell, hash_prefix)) {
    return(sub(hash_prefix, "", cell, fixed = TRUE))
  }
  if (endsWith(cell, suffix)) {
    return(substr(cell, 1L, nchar(cell) - nchar(suffix)))
  }
  cell
}

archr_fragment_barcode <- function(cell, sample, suffix = NULL) {
  barcode <- archr_strip_sample_from_cell(cell, sample)
  if (!is.null(suffix) && length(suffix) == 1L && nzchar(suffix)) {
    # Some multi-library exports encode the library in the barcode suffix
    # (for example ...-2-s1d2), while the corresponding fragment file uses
    # the original 10x barcode suffix (...-1). Normalize only when the
    # dataset explicitly declares this rule.
    barcode <- sub("-[0-9]+$", "", barcode)
    barcode <- paste0(barcode, suffix)
  }
  barcode
}

archr_arrow_fragment_chromosomes <- function(arrow_files) {
  chromosome_sets <- lapply(arrow_files, function(arrow_file) {
    h5 <- rhdf5::h5ls(arrow_file, recursive = TRUE)
    as.character(h5$name[
      h5$group == "/Fragments" & h5$otype == "H5I_GROUP"
    ])
  })
  chromosome_sets <- chromosome_sets[lengths(chromosome_sets) > 0L]
  if (!length(chromosome_sets)) {
    stop("Could not find any /Fragments chromosomes in ArchR Arrow files",
         call. = FALSE)
  }
  Reduce(intersect, chromosome_sets)
}

archr_real_make_sce <- function(rna, gene_annotation, sample_name,
                                fragment_barcode_suffix = NULL) {
  sample_name <- as.character(sample_name)
  if (length(sample_name) == 1L) {
    sample_name <- rep(sample_name, ncol(rna))
  }
  if (length(sample_name) != ncol(rna)) {
    stop("ArchR sample-name mapping does not match RNA cells", call. = FALSE)
  }
  genes <- gene_annotation$genes
  gene_metadata <- as.data.frame(S4Vectors::mcols(genes))
  symbol_column <- intersect(c("symbol", "name", "gene_name"), names(gene_metadata))
  if (!length(symbol_column)) {
    stop("ArchR gene annotation has no gene-symbol column", call. = FALSE)
  }
  gene_symbols <- as.character(gene_metadata[[symbol_column[[1L]]]])
  gene_index <- match(rownames(rna), gene_symbols)
  keep <- !is.na(gene_index)
  if (!any(keep)) {
    stop("No RNA features matched the ArchR gene annotation", call. = FALSE)
  }
  rna <- rna[keep, , drop = FALSE]
  matched_genes <- genes[gene_index[keep]]
  S4Vectors::mcols(matched_genes)$name <- rownames(rna)
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = Matrix::Matrix(rna, sparse = TRUE)),
    rowRanges = matched_genes
  )
  rownames(sce) <- rownames(rna)
  # Arrow files use sample#barcode cell names.  The input RNA/ATAC matrices
  # use sample_barcode names for multi-sample datasets, so remove the sample
  # prefix before constructing the SCE name.  Keeping sample_barcode here
  # would produce sample#sample_barcode and make addGeneExpressionMatrix()
  # report zero overlapping cells.
  raw_barcodes <- vapply(seq_along(sample_name), function(i) {
    archr_fragment_barcode(
      colnames(rna)[[i]], sample_name[[i]], fragment_barcode_suffix
    )
  }, character(1))
  colnames(sce) <- paste0(sample_name, "#", raw_barcodes)
  SummarizedExperiment::rowData(sce)$name <- rownames(rna)
  sce
}

archr_real_link_results <- function(project, pairs, fdr) {
  links <- ArchR::getPeak2GeneLinks(
    project, corCutOff = 0, FDRCutOff = 1,
    # Current ArchR versions validate these as numeric inputs. Zero keeps
    # the result extraction from adding an extra variance filter; FDR is
    # applied below so that all methods use the same final cutoff.
    varCutOffATAC = 0, varCutOffRNA = 0, returnLoops = FALSE
  )
  results <- data.frame(
    gene = pairs$gene, peak = pairs$peak,
    Correlation = NA_real_, FDR = NA_real_, pval = NA_real_,
    padj = NA_real_, stringsAsFactors = FALSE
  )
  if (is.null(links) || !nrow(links)) {
    return(list(results = results, significant = results[FALSE, , drop = FALSE]))
  }

  # idxATAC/idxRNA refer to the feature sets stored in the links object,
  # rather than to the full project gene annotation.  Using getGenes(project)
  # here can silently produce a different ordering and zero pair matches.
  link_peak_set <- S4Vectors::metadata(links)$peakSet
  link_gene_set <- S4Vectors::metadata(links)$geneSet
  if (is.null(link_peak_set) || is.null(link_gene_set)) {
    stop("ArchR link output is missing its peakSet or geneSet metadata", call. = FALSE)
  }
  peak_names <- names(link_peak_set)
  if (!length(peak_names)) {
    peak_names <- paste(
      as.character(GenomicRanges::seqnames(link_peak_set)),
      GenomicRanges::start(link_peak_set),
      GenomicRanges::end(link_peak_set), sep = "-"
    )
  }
  peak_metadata <- as.data.frame(S4Vectors::mcols(link_peak_set))
  if ("name" %in% names(peak_metadata)) {
    stored_names <- as.character(peak_metadata$name)
    if (all(!is.na(stored_names)) && all(nzchar(stored_names))) {
      peak_names <- stored_names
    }
  }
  gene_metadata <- as.data.frame(S4Vectors::mcols(link_gene_set))
  symbol_column <- intersect(c("name", "symbol", "gene_name"), names(gene_metadata))
  gene_names <- if (length(symbol_column)) {
    as.character(gene_metadata[[symbol_column[[1L]]]])
  } else {
    names(link_gene_set)
  }
  if (is.null(peak_names) || is.null(gene_names)) {
    stop("ArchR link output is missing peak or gene names", call. = FALSE)
  }

  linked_peak <- peak_names[as.integer(links$idxATAC)]
  linked_gene <- gene_names[as.integer(links$idxRNA)]
  link_key <- paste(linked_gene, linked_peak, sep = "\r")
  pair_key <- paste(pairs$gene, pairs$peak, sep = "\r")
  result_index <- match(pair_key, link_key)
  results$Correlation <- as.numeric(links$Correlation[result_index])
  results$FDR <- as.numeric(links$FDR[result_index])
  if ("Pval" %in% names(links)) {
    results$pval <- as.numeric(links$Pval[result_index])
  }
  results$padj <- results$FDR
  significant <- results[
    is.finite(results$padj) & results$padj <= fdr,
    , drop = FALSE
  ]
  list(results = results, significant = significant)
}

archr_real_annotations <- function(genome_name) {
  ArchR::addArchRGenome(genome_name, install = FALSE)
  default_annotations <- tryCatch(
    list(
      gene = ArchR::getGeneAnnotation(),
      genome = ArchR::getGenomeAnnotation()
    ),
    error = function(e) NULL
  )
  if (!is.null(default_annotations)) return(default_annotations)

  # Some ArchR installations have the package data index but are missing the
  # serialized geneAnnoHg38 object.  Rebuild the same type of annotation from
  # the Bioconductor EnsDb package instead of relying on ArchR's data cache.
  ensdb_package <- switch(
    genome_name,
    hg38 = "EnsDb.Hsapiens.v86",
    hg19 = "EnsDb.Hsapiens.v75",
    NULL
  )
  if (is.null(ensdb_package) ||
      !requireNamespace(ensdb_package, quietly = TRUE) ||
      !requireNamespace("ensembldb", quietly = TRUE)) {
    stop(
      "ArchR built-in gene annotation is unavailable for ", genome_name,
      "; install the matching EnsDb and ensembldb packages",
      call. = FALSE
    )
  }
  ensdb <- getExportedValue(ensdb_package, ensdb_package)
  genes <- ensembldb::genes(ensdb, return.type = "GRanges")
  exons <- ensembldb::exons(ensdb, return.type = "GRanges")
  transcripts <- ensembldb::transcripts(ensdb, return.type = "GRanges")
  to_ucsc_seqnames <- function(gr) {
    old <- GenomeInfoDb::seqlevels(gr)
    new <- old
    primary <- old %in% c(as.character(seq_len(22L)), "X", "Y")
    new[primary] <- paste0("chr", old[primary])
    new[old %in% c("M", "MT")] <- "chrM"
    GenomeInfoDb::renameSeqlevels(gr, setNames(new, old))
  }
  genes <- to_ucsc_seqnames(genes)
  exons <- to_ucsc_seqnames(exons)
  transcripts <- to_ucsc_seqnames(transcripts)
  gene_metadata <- as.data.frame(S4Vectors::mcols(genes))
  gene_id <- if ("gene_id" %in% names(gene_metadata)) {
    as.character(gene_metadata$gene_id)
  } else {
    names(genes)
  }
  gene_symbol <- if ("gene_name" %in% names(gene_metadata)) {
    as.character(gene_metadata$gene_name)
  } else {
    gene_id
  }
  missing_symbol <- is.na(gene_symbol) | !nzchar(gene_symbol)
  gene_symbol[missing_symbol] <- gene_id[missing_symbol]
  S4Vectors::mcols(genes)$symbol <- gene_symbol
  names(genes) <- gene_symbol
  tss <- unique(GenomicRanges::resize(transcripts, width = 1L, fix = "start"))
  genome <- ArchR::createGenomeAnnotation(
    genome = genome_name, blacklist = GenomicRanges::GRanges(),
    filter = TRUE, filterChr = c("chrM")
  )
  list(
    gene = ArchR::createGeneAnnotation(
      genes = genes, exons = exons, TSS = tss
    ),
    genome = genome
  )
}

#' Run ArchR Peak2GeneLinks on one real multiome cell type.
run_archr_real <- function(rna, atac, pairs, keep_cells, fragment_file = NULL,
                           arrow_files = NULL, sample_name, genome_name,
                           output_dir, fdr, n_cores = 1L, seed = 1L,
                           max_dist = 250000L, k = 100L,
                           fragment_barcode_suffix = NULL) {
  archr_requirements()
  use_existing_arrows <- !is.null(arrow_files)
  fragment_files <- if (is.null(fragment_file)) {
    character()
  } else {
    as.character(fragment_file)
  }
  arrow_files <- if (is.null(arrow_files)) {
    character()
  } else {
    as.character(arrow_files)
  }
  if (use_existing_arrows && length(fragment_files)) {
    stop("ArchR cannot use fragments and existing Arrow files together",
         call. = FALSE)
  }
  if (!use_existing_arrows && !length(fragment_files)) {
    stop("ArchR requires fragment files or existing Arrow files", call. = FALSE)
  }
  sample_names <- as.character(sample_name)
  input_files <- if (use_existing_arrows) arrow_files else fragment_files
  if (!length(input_files) || any(!nzchar(input_files))) {
    stop("ArchR input files cannot be empty", call. = FALSE)
  }
  if (length(sample_names) == 1L && length(input_files) > 1L) {
    sample_names <- sub(
      "(_atac_fragments[.]tsv[.]gz|[.]arrow)$", "", basename(input_files)
    )
  }
  if (length(sample_names) != length(input_files) ||
      any(!nzchar(sample_names)) || anyDuplicated(sample_names)) {
    stop(
      "ArchR input files and sample names must have the same unique length",
      call. = FALSE
    )
  }
  if (any(!file.exists(input_files))) {
    stop(
      "ArchR input file not found: ",
      paste(input_files[!file.exists(input_files)], collapse = ", "),
      call. = FALSE
    )
  }
  if (!use_existing_arrows) {
    fragment_indexes <- paste0(fragment_files, ".tbi")
    if (any(!file.exists(fragment_indexes))) {
      stop(
        "ArchR fragment index not found: ",
        paste(fragment_indexes[!file.exists(fragment_indexes)], collapse = ", "),
        call. = FALSE
      )
    }
  }
  # Fragment files commonly store bare 10x barcodes (e.g. AAAC...-1),
  # whereas some matrix exports store sample_barcode names.  Support both
  # forms.  ArchR adds sample# internally when it creates the Arrow file.
  selected_counts <- vapply(
    sample_names,
    function(sample) {
      sum(vapply(
        keep_cells, archr_cell_matches_sample, logical(1), sample = sample
      ))
    },
    integer(1)
  )
  if (!any(selected_counts > 0L) && length(fragment_files) == 1L) {
    selected_counts[[1L]] <- length(keep_cells)
  }
  if (!any(selected_counts > 0L)) {
    stop("No selected cells match any ArchR fragment sample", call. = FALSE)
  }
  if (any(selected_counts == 0L)) {
    keep_samples <- selected_counts > 0L
    input_files <- input_files[keep_samples]
    if (use_existing_arrows) {
      arrow_files <- input_files
    } else {
      fragment_files <- input_files
    }
    sample_names <- sample_names[keep_samples]
  }
  if (!genome_name %in% c("hg19", "hg38")) {
    stop("ArchR real-data support currently expects hg19 or hg38", call. = FALSE)
  }
  if (length(keep_cells) < 4L) {
    stop("ArchR requires at least four selected cells", call. = FALSE)
  }

  work_directory <- file.path(output_dir, "archr_work")
  if (grepl("[[:space:]]", work_directory)) {
    # ArchRProject rejects outputDirectory paths containing spaces.  Keep the
    # benchmark results in their human-readable output directory, but place
    # ArchR's internal Arrow/project workspace in a sibling path without
    # spaces.
    safe_label <- gsub(
      "[^A-Za-z0-9_.-]+", "_", basename(dirname(output_dir))
    )
    safe_method <- gsub("[^A-Za-z0-9_.-]+", "_", basename(output_dir))
    work_directory <- file.path(
      dirname(dirname(output_dir)),
      paste0(safe_label, "_", safe_method, "_work")
    )
  }
  arrow_directory <- file.path(work_directory, "arrow")
  project_directory <- file.path(work_directory, "project")
  dir.create(arrow_directory, recursive = TRUE, showWarnings = FALSE)
  dir.create(project_directory, recursive = TRUE, showWarnings = FALSE)

  ArchR::addArchRThreads(threads = as.integer(n_cores))
  ArchR::addArchRLocking(locking = TRUE)
  annotations <- archr_real_annotations(genome_name)
  gene_annotation <- annotations$gene
  genome_annotation <- annotations$genome
  peak_ranges <- archr_peak_ranges(list(atac = atac))
  cell_sample_names <- vapply(
    colnames(rna),
    function(cell) {
      hits <- sample_names[vapply(
        sample_names, archr_cell_matches_sample, logical(1), cell = cell
      )]
      if (!length(hits) && length(sample_names) == 1L) {
        hits <- sample_names
      }
      if (length(hits) != 1L) {
        stop(
          "Could not map RNA cell to exactly one ArchR fragment sample: ",
          cell, call. = FALSE
        )
      }
      hits[[1L]]
    },
    character(1)
  )
  se_rna <- archr_real_make_sce(
    rna, gene_annotation, cell_sample_names, fragment_barcode_suffix
  )

  valid_barcodes <- lapply(sample_names, function(sample) {
    selected <- keep_cells[vapply(
      keep_cells, archr_cell_matches_sample, logical(1), sample = sample
    )]
    if (!length(selected) && length(fragment_files) == 1L) {
      # Some single-library datasets store raw 10x barcodes without a sample
      # prefix or suffix (for example Pituitary).
      selected <- keep_cells
    }
    if (!length(selected)) {
      stop(
        "No selected cells match ArchR fragment sample ", sample,
        call. = FALSE
      )
    }
    vapply(
      selected, archr_fragment_barcode, character(1), sample = sample,
      suffix = fragment_barcode_suffix
    )
  })
  names(valid_barcodes) <- sample_names
  if (length(fragment_files) == 1L) {
    valid_barcodes <- valid_barcodes[[1L]]
  }

  if (!use_existing_arrows) {
    arrow_files <- ArchR::createArrowFiles(
      inputFiles = fragment_files, sampleNames = sample_names,
      outputNames = file.path(arrow_directory, sample_names),
      # For multiple inputs ArchR 1.0.3 expects a sample-name keyed list;
      # for one input it accepts the barcode vector directly.
      validBarcodes = valid_barcodes,
      geneAnnotation = gene_annotation, genomeAnnotation = genome_annotation,
      minTSS = 0, minFrags = 0, maxFrags = 1e9,
      addTileMat = FALSE, addGeneScoreMat = FALSE,
      filterFrags = 0, filterTSS = 0,
      threads = as.integer(n_cores), subThreading = FALSE,
      force = TRUE, verbose = FALSE, cleanTmp = FALSE,
      logFile = file.path(work_directory, "createArrowFiles.log")
    )
    # ArchR 1.0.3 can return a list-like object here when only one Arrow file
    # is created.  ArchRProject later passes sampleColData$ArrowFiles directly
    # to file.exists(), which requires an ordinary character vector.
    arrow_files <- as.character(unlist(arrow_files, use.names = FALSE))
    arrow_files <- arrow_files[nzchar(arrow_files)]
    if (!length(arrow_files) || any(!file.exists(arrow_files))) {
      stop(
        "ArchR did not create a valid Arrow file: ",
        paste(arrow_files[!file.exists(arrow_files)], collapse = ", "),
        call. = FALSE
      )
    }
  } else {
    # Existing Arrow files already contain the fragment index and PeakMatrix.
    # Keep their paths unchanged; ArchRProject will attach the configured
    # sample names to the cell names (sample#barcode).
    arrow_files <- as.character(arrow_files)
  }
  if (!use_existing_arrows) {
    arrow_chromosomes <- archr_arrow_fragment_chromosomes(arrow_files)
    peak_chromosomes <- as.character(GenomicRanges::seqnames(peak_ranges))
    keep_peaks <- peak_chromosomes %in% arrow_chromosomes
    if (!all(keep_peaks)) {
      message(
        "ArchR: dropping ", sum(!keep_peaks),
        " peaks on chromosomes absent from Arrow fragment files: ",
        paste(sort(unique(peak_chromosomes[!keep_peaks])), collapse = ",")
      )
      peak_ranges <- peak_ranges[keep_peaks]
    }
    if (!length(peak_ranges)) {
      stop("No ATAC peaks remain on chromosomes present in Arrow fragments",
           call. = FALSE)
    }
  }
  project <- ArchR::ArchRProject(
    ArrowFiles = arrow_files, outputDirectory = project_directory,
    copyArrows = FALSE, geneAnnotation = gene_annotation,
    genomeAnnotation = genome_annotation, showLogo = FALSE
  )
  if (use_existing_arrows) {
    # Arrow files retain the PeakMatrix but ArchRProject() does not always
    # reconstruct the project-level peakSet metadata from that matrix.
    # Attach the repository's exact peak coordinates before subsetting and
    # running Peak2GeneLinks.
    project <- ArchR::addPeakSet(project, peakSet = peak_ranges, force = TRUE)
    project_cells <- intersect(keep_cells, ArchR::getCellNames(project))
    if (length(project_cells) < 4L) {
      stop(
        "Existing Arrow files contain fewer than four selected cells; mapped ",
        length(project_cells), call. = FALSE
      )
    }
    project <- ArchR::subsetArchRProject(
      ArchRProj = project, cells = project_cells,
      outputDirectory = file.path(work_directory, "project_subset"),
      dropCells = TRUE, force = TRUE, threads = as.integer(n_cores),
      logFile = file.path(work_directory, "subsetArchRProject.log")
    )
    keep_cells <- project_cells
    se_rna <- se_rna[, intersect(colnames(se_rna), project_cells), drop = FALSE]
    if (!length(colnames(se_rna)) ||
        !all(project_cells %in% colnames(se_rna))) {
      stop("RNA cells do not fully overlap the existing ArchR Arrow cells",
           call. = FALSE)
    }
  } else {
    project <- ArchR::addPeakSet(project, peakSet = peak_ranges, force = TRUE)
    project <- ArchR::addPeakMatrix(
      project, force = TRUE, threads = as.integer(n_cores),
      verbose = FALSE, logFile = file.path(work_directory, "addPeakMatrix.log")
    )
  }
  project <- ArchR::addGeneExpressionMatrix(
    input = project, seRNA = se_rna, strictMatch = TRUE,
    scaleTo = 10000, force = TRUE, threads = as.integer(n_cores),
    verbose = FALSE,
    logFile = file.path(work_directory, "addGeneExpressionMatrix.log")
  )

  n_cells <- length(keep_cells)
  n_features <- nrow(atac)
  dims <- archr_lsi_dims(n_cells, n_features)
  project <- ArchR::addIterativeLSI(
    ArchRProj = project, useMatrix = "PeakMatrix",
    name = "IterativeLSI", iterations = 1L,
    varFeatures = archr_var_features(n_features),
    dimsToUse = seq_len(dims), outlierQuantiles = NULL,
    filterBias = FALSE, sampleCellsPre = NULL, sampleCellsFinal = NULL,
    keep0lsi = TRUE,
    clusterParams = list(
      resolution = 0.2, sampleCells = n_cells, n.start = 1L
    ),
    threads = as.integer(n_cores), seed = as.integer(seed), force = TRUE,
    verbose = FALSE, logFile = file.path(work_directory, "addIterativeLSI.log")
  )
  group_k <- as.integer(max(2L, min(k, floor(n_cells / 4L), n_cells - 1L)))
  project <- ArchR::addPeak2GeneLinks(
    ArchRProj = project, reducedDims = "IterativeLSI",
    useMatrix = "GeneExpressionMatrix", dimsToUse = seq_len(dims),
    k = group_k, knnIteration = min(500L, n_cells), overlapCutoff = 0.8,
    maxDist = as.integer(max_dist), scaleTo = 10000, log2Norm = TRUE,
    predictionCutoff = 0.0, addEmpiricalPval = FALSE,
    addPermutedPval = FALSE, seed = as.integer(seed),
    threads = as.integer(n_cores), verbose = FALSE,
    logFile = file.path(work_directory, "addPeak2GeneLinks.log")
  )

  links <- archr_real_link_results(project, pairs, fdr)
  list(project = project, results = links$results,
       significant = links$significant)
}

archr_peak_ranges <- function(simulation) {
  peaks <- rownames(simulation$atac)
  if (is.null(peaks) && methods::is(simulation$atac, "Matrix")) {
    peaks <- methods::slot(simulation$atac, "Dimnames")[[1L]]
  }
  if (is.null(peaks)) peaks <- dimnames(simulation$atac)[[1L]]
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
    varCutOffATAC = 0, varCutOffRNA = 0,
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
    work_directory <- tempfile("calibre-archr-")
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
    force = TRUE, verbose = FALSE, cleanTmp = FALSE,
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
  dims <- archr_lsi_dims(simulation$n_cells, simulation$n_pairs, cap = 10L)
  project <- ArchR::addIterativeLSI(
    ArchRProj = project, useMatrix = "PeakMatrix",
    name = "IterativeLSI", iterations = 1L,
    varFeatures = archr_var_features(simulation$n_pairs),
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
