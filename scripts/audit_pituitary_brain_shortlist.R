#!/usr/bin/env Rscript

# Exact-pair, Poisson-only audit of the pituitary/brain shortlist.  Brain2 ALL
# uses the pooled EpiLink model adjusted for condition and compares it with the
# pooled scMultiMap and Signac outputs from the same cells.

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[[1L]] else ".", mustWork = TRUE)
out <- file.path(root, "results", "validation", "functional", "GWAS_loci",
                 "pituitary_brain_shortlist_audit.tsv")

cfg <- data.table::data.table(
  candidate = c("FSHB_rs11031006", "GH1_promoter_rs11568828", "GH1_coding_rs5388",
                "ESR1_endometriosis", "POU1F1_height", "LHB_testosterone",
                "CLU_rs1532278", "SORL1_rs2276412", "GRIK5_AD", "ZEB2_migraine",
                "APOE_rs7256200", "TSPAN14_rs7922621"),
  dataset = c(rep("Pituitary", 6), rep("Brain2", 6)),
  celltype = c("Gonadotropes", "Somatotropes", "Somatotropes", "Gonadotropes",
               "Lactotropes", "Gonadotropes", "Astrocytes", "Microglia",
               "Astrocytes", "Astrocytes", "Astrocytes", "Microglia"),
  state = c(rep("standard", 6), rep("ALL", 6)),
  gene = c("FSHB", "GH1", "GH1", "ESR1", "POU1F1", "LHB", "CLU", "SORL1",
           "GRIK5", "ZEB2", "APOE", "TSPAN14"),
  peak = c("chr11-30203591-30205849", "chr17-63918218-63919459",
           "chr17-63915547-63917975", "chr6-151493400-151494947",
           "chr3-87276200-87276766", "chr19-49512181-49514450",
           "chr8-27608401-27609282", "chr11-121589776-121590720",
           "chr19-42069748-42070671", "chr2-144510178-144511104",
           "chr19-44912591-44913496", "chr10-80491191-80492046"),
  variant = c("rs11031006", "rs11568828", "rs5388", "6:151494876:T:C",
              "3:87276603", "rs111981233", "rs1532278", "rs2276412",
              "19:42069985:C:T", "2:144510687:CCTCT:C", "rs7256200", "rs7922621"),
  celltype_evidence = c("matched", "matched", "matched", "matched", "POU1F1-lineage",
                        "matched", "conflicts with exact-SNP experiment", "matched",
                        "plausible but unvalidated", "uncertain", "matched", "matched"),
  external_evidence = c(
    "allele-specific SF1 binding/reporter and knock-in mouse",
    "promoter haplotype reporter/EMSA, but exact Somatotrope link is null",
    "coding missense variant; clinical/GH biology, not enhancer-to-gene validation",
    "endometriosis fine-mapping/eQTL/HiChIP near ESR1; no pituitary allele perturbation",
    "POU1F1 pituitary lineage and height biology; no exact-allele validation",
    "correct gonadotrope hormone biology; variant region also supports FCGRT and lacks LHB perturbation",
    "exact isogenic CRISPR supports CLU in excitatory neurons and is negative in astrocytes",
    "Yang microglia pcHi-C target plus SORL1 microglial KO biology; rs2276412 is synonymous coding",
    "very high FinnGen AD PIP but no exact allele-to-GRIK5 experiment found",
    "reproducible EpiLink link in astrocytes but no exact allele-to-ZEB2 validation found",
    "ME1 haplotype reporter evidence; older study found no effect in U118 astrocytes",
    "direct microglial prime-edit/CRISPR evidence, but EpiLink pair is not significant")
)

paths <- function(ds, ct, state) {
  base <- file.path(root, "results", "real", ds, ct)
  if (ds == "Brain2") {
    base <- file.path(base, state)
    ep_dir <- if (ct == "Microglia") "epilink_poisson_condition_effect" else
      "epilink_poisson_condition_effect"
  } else ep_dir <- "epilink_poisson"
  list(ep = file.path(base, ep_dir, "all_pairs.tsv.gz"),
       sm = file.path(base, "scmultimap", "all_pairs.tsv.gz"),
       sg = file.path(base, "signac", "all_pairs.tsv.gz"))
}

one <- function(row) {
  p <- paths(row$dataset, row$celltype, row$state)
  getrow <- function(path, cols) {
    if (!file.exists(path)) return(data.table())
    h <- names(fread(path, nrows = 0L, showProgress = FALSE))
    # Stream only lines containing the exact peak.  Several all-pairs files are
    # large enough that repeatedly loading them in full is wasteful.
    cmd <- paste("zgrep -F", shQuote(row$peak), shQuote(path))
    x <- tryCatch(fread(cmd = cmd, header = FALSE, col.names = h,
                        showProgress = FALSE), error = function(e) data.table())
    if (!nrow(x)) return(data.table())
    x <- x[, intersect(cols, h), with = FALSE]
    x[gene == row$gene & peak == row$peak]
  }
  ep <- getrow(p$ep, c("gene", "peak", "beta", "pval", "pvalue", "padj"))
  sm <- getrow(p$sm, c("gene", "peak", "pval", "padj", "cor"))
  sg <- getrow(p$sg, c("gene", "peak", "pvalue", "padj", "pvalue_uncorrected",
                        "padj_uncorrected", "score"))
  data.table(
    epilink_beta = if (nrow(ep) && "beta" %in% names(ep)) ep$beta[1] else NA_real_,
    epilink_p = if (nrow(ep)) ep[[intersect(c("pval", "pvalue"), names(ep))[1]]][1] else NA_real_,
    epilink_fdr = if (nrow(ep)) ep$padj[1] else NA_real_,
    scmultimap_fdr = if (nrow(sm)) sm$padj[1] else NA_real_,
    signac_corrected_fdr = if (nrow(sg)) sg$padj[1] else NA_real_,
    signac_uncorrected_fdr = if (nrow(sg)) sg$padj_uncorrected[1] else NA_real_
  )
}

stats <- rbindlist(lapply(seq_len(nrow(cfg)), function(i) one(cfg[i])))
ans <- cbind(cfg, stats)
ans[, epilink_only := is.finite(epilink_fdr) & epilink_fdr < 0.1 &
       (!is.finite(scmultimap_fdr) | scmultimap_fdr >= 0.1) &
       (!is.finite(signac_corrected_fdr) | signac_corrected_fdr >= 0.1) &
       (!is.finite(signac_uncorrected_fdr) | signac_uncorrected_fdr >= 0.1)]
ans[, verdict := fifelse(candidate == "FSHB_rs11031006", "positive-control; strongest pituitary",
                  fifelse(candidate == "CLU_rs1532278", "statistical positive; cell-type contradiction",
                  fifelse(candidate == "SORL1_rs2276412", "orthogonal-support candidate; not enhancer-clean",
                  fifelse(candidate == "ZEB2_migraine" & epilink_only, "novel hypothesis; needs experiment",
                  fifelse(epilink_only, "hypothesis; evidence incomplete", "reject: EpiLink criterion fails")))))]
fwrite(ans, out, sep = "\t")
message("Wrote ", out)
