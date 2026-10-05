#!/usr/bin/env Rscript

# Standardise FinnGen R12, RA, T1D and IBD fine-mapping data for the
# SCENT-style causal-variant enrichment benchmark.

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(GenomicRanges)
  library(IRanges)
  library(rtracklayer)
})

args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[[1]] else ".", mustWork = TRUE)
gwas <- file.path(root, "results", "validation", "GWAS")
ref_root <- file.path(gwas, "reference")
chain_file <- file.path(ref_root, "liftover", "hg19ToHg38.over.chain.gz")

read_chain_gz <- function(path) {
  tmp <- tempfile(fileext = ".chain")
  writeLines(readLines(gzfile(path), warn = FALSE), tmp)
  import.chain(tmp)
}
chain19to38 <- read_chain_gz(chain_file)

empty_variants <- function() data.table(
  source=character(), trait=character(), locus=character(), chr=character(),
  pos=integer(), pip=numeric(), consequence=character(), variant_id=character())

lift_points <- function(x) {
  if (!nrow(x)) return(x[, .(source, trait, locus, chr, pos, pip, consequence, variant_id)])
  gr <- GRanges(x$chr, IRanges(x$pos, x$pos))
  lo <- liftOver(gr, chain19to38)
  ok <- lengths(lo) == 1L
  ans <- copy(x[ok])
  one <- unlist(lo[ok], use.names = FALSE)
  ans[, `:=`(chr = as.character(seqnames(one)), pos = as.integer(start(one)))]
  attr(ans, "lift_n_input") <- nrow(x)
  attr(ans, "lift_n_output") <- nrow(ans)
  ans
}

coding_re <- paste(c(
  "missense", "nonsynonymous", "stopgain", "stop_gained", "stoploss",
  "stop_lost", "start_lost", "frameshift", "splice", "splicing"
), collapse = "|")

make_loci_from_variants <- function(v, flank = 500000L) {
  # External papers do not all publish the exact analysis-window boundaries.
  # For those sources use a reproducible +/-500 kb window around the
  # highest-PIP variant; IBD supplies explicit region boundaries separately.
  v[, {
    lead <- which.max(pip)
    list(chr=chr[lead], start=max(1L, pos[lead]-flank), end=pos[lead]+flank)
  }, by=.(source, trait, locus)]
}

cluster_signal_loci <- function(v, flank=500000L) {
  leads <- v[, {i <- which.max(pip); list(chr=chr[i], lead_pos=pos[i])},
             by=.(source,trait,old_locus=locus)]
  win <- GRanges(leads$chr, IRanges(pmax(1L,leads$lead_pos-flank), leads$lead_pos+flank))
  red <- reduce(win, ignore.strand=TRUE)
  hit <- findOverlaps(win, red, select="first", ignore.strand=TRUE)
  map <- leads[, .(source,trait,old_locus)]
  map[, locus := paste0(source,"::locus_",sprintf("%03d",hit))]
  v2 <- copy(v)
  setnames(v2, "locus", "old_locus")
  v2 <- merge(v2, map, by=c("source","trait","old_locus"), all.x=TRUE)
  v2[, old_locus := NULL]
  loci <- data.table(source=unique(v$source), trait=unique(v$trait),
                     locus=paste0(unique(v$source),"::locus_",sprintf("%03d",seq_along(red))),
                     chr=as.character(seqnames(red)), start=start(red), end=end(red))
  list(variants=v2, loci=loci)
}

apply_qc <- function(v, loci, consequence_available = TRUE) {
  v <- unique(v[!is.na(chr) & !is.na(pos) & !is.na(pip)], by=c("source","trait","locus","chr","pos"))
  v <- v[!(chr == "chr6" & pos >= 25000000L & pos <= 36000000L)]
  loci <- loci[!(chr == "chr6" & end >= 25000000L & start <= 36000000L)]
  bad <- character()
  if (consequence_available) {
    bad <- unique(v[pip > 0.5 & grepl(coding_re, consequence, ignore.case=TRUE), locus])
    v <- v[!locus %chin% bad]
    loci <- loci[!locus %chin% bad]
  }
  causal <- v[pip > 0.2]
  list(variants=v, causal=causal, loci=unique(loci), removed_loci=bad,
       consequence_qc=consequence_available)
}

write_reference <- function(source, obj, raw_description, build_in, locus_definition,
                            lift_input=NA_integer_, lift_output=NA_integer_) {
  if (is.null(lift_input) || !length(lift_input)) lift_input <- NA_integer_
  if (is.null(lift_output) || !length(lift_output)) lift_output <- NA_integer_
  out <- file.path(ref_root, source, "processed")
  dir.create(out, recursive=TRUE, showWarnings=FALSE)
  fwrite(obj$variants, file.path(out, "variants_all.tsv.gz"), sep="\t")
  fwrite(obj$causal, file.path(out, "causal_variants_PIPgt0.2.tsv.gz"), sep="\t")
  fwrite(obj$loci, file.path(out, "analysis_loci.hg38.tsv.gz"), sep="\t")
  manifest <- data.table(
    source=source, raw_description=raw_description, input_build=build_in,
    analysis_build="GRCh38", pip_operator=">", pip_threshold=0.2,
    coding_splicing_locus_qc=obj$consequence_qc,
    coding_splicing_qc_rule="remove locus when coding/splicing variant PIP > 0.5",
    mhc_excluded="chr6:25-36Mb", locus_definition=locus_definition,
    n_variants=nrow(obj$variants), n_causal=nrow(obj$causal),
    n_traits=uniqueN(obj$causal$trait), n_loci=nrow(obj$loci),
    n_loci_removed_consequence=length(obj$removed_loci),
    lift_input=lift_input, lift_output=lift_output,
    lift_fraction=if (!is.na(lift_input) && lift_input > 0) lift_output/lift_input else NA_real_
  )
  fwrite(manifest, file.path(out, "reference_manifest.tsv"), sep="\t")
  saveRDS(list(source=source, causal=obj$causal, loci=obj$loci, manifest=manifest),
          file.path(out, "scent_reference.rds"), compress="xz")
}

message("Preparing FinnGen R12 (>1,000 cases)")
fg_dir <- file.path(ref_root, "FinnGen")
fg_manifest <- fread(file.path(fg_dir, "finngen_R12_manifest.tsv"))
eligible <- fg_manifest[num_cases > 1000]
fg_files <- file.path(fg_dir, "raw", "summary",
                      paste0("finngen_R12_", eligible$phenocode, ".SUSIE.snp.filter.tsv"))
fg_files <- fg_files[file.exists(fg_files)]
fg_list <- vector("list", length(fg_files))
for (i in seq_along(fg_files)) {
  x <- fread(fg_files[[i]], showProgress=FALSE)
  if (!nrow(x)) next
  x <- x[!is.na(position) & !is.na(cs_specific_prob)]
  if (!nrow(x)) next
  x[, pip := as.numeric(cs_specific_prob)]
  x[, `:=`(
    source="FinnGen", trait=as.character(trait), locus=paste(trait, region, sep="::"),
    chr=ifelse(grepl("^chr", chromosome), chromosome, paste0("chr", chromosome)),
    pos=as.integer(position), consequence=as.character(most_severe), variant_id=as.character(v)
  )]
  fg_list[[i]] <- x[, .(source,trait,locus,chr,pos,pip,consequence,variant_id)]
}
fg <- rbindlist(fg_list, fill=TRUE)
fg <- fg[, .(pip=max(pip), consequence=consequence[which.max(pip)],
             variant_id=variant_id[which.max(pip)]),
         by=.(source,trait,locus,chr,pos)]
parse_region <- function(z) {
  m <- regexec("^(chr[^:]+):([0-9]+)-([0-9]+)$", z)
  p <- regmatches(z,m)
  data.table(chr=vapply(p, `[`, character(1), 2),
             start=as.integer(vapply(p, `[`, character(1), 3)),
             end=as.integer(vapply(p, `[`, character(1), 4)))
}
fg_loci <- unique(fg[, .(source,trait,locus)])
fg_loci <- cbind(fg_loci, parse_region(sub("^.*::", "", fg_loci$locus)))
fg_obj <- apply_qc(fg, fg_loci, TRUE)
write_reference("FinnGen", fg_obj,
  "FinnGen release 12 public SuSiE 95% credible-set files; endpoints with >1,000 cases",
  "GRCh38", "published FinnGen fine-mapping region")

message("Preparing rheumatoid arthritis")
ra_file <- file.path(ref_root,"RA","Ishigaki_2022_supplementary_tables.xlsx")
ra0 <- as.data.table(read_excel(ra_file, sheet="ST-6", skip=5, col_names=FALSE))
ra_vid <- tstrsplit(as.character(ra0$...7), ":", fixed=TRUE)
ra <- ra0[, .(
  source="RA", trait="Rheumatoid arthritis",
  locus=paste0("RA::", as.character(...3)),
  chr=paste0("chr", ra_vid[[1]]), pos=as.integer(ra_vid[[2]]),
  pip=as.numeric(...13),
  consequence=ifelse(!is.na(...17) & as.character(...17) != ".",
                     paste("nonsynonymous",...15,...17,sep=";"),
                     paste(...15,...17,sep=";")),
  variant_id=as.character(...7)
)]
ra <- ra[!is.na(pos) & !is.na(pip)]
ra_lift <- lift_points(ra)
ra_lift_n_input <- attr(ra_lift,"lift_n_input")
ra_lift_n_output <- attr(ra_lift,"lift_n_output")
ra_loci <- make_loci_from_variants(ra_lift)
ra_obj <- apply_qc(ra_lift, ra_loci, TRUE)
write_reference("RA", ra_obj,
  "Ishigaki et al. 2022 Supplementary Table 6 multi-ancestry fine-mapping",
  "GRCh37/hg19", "+/-500 kb around highest-PIP variant per published lead-association locus",
  ra_lift_n_input, ra_lift_n_output)

message("Preparing type 1 diabetes")
t1_file <- file.path(ref_root,"T1D","Chiou_2021_T1D_99pct_credible_sets.xlsx")
t10 <- as.data.table(read_excel(t1_file, sheet="Supplementary Data 1", skip=2))
t1 <- t10[, .(
  source="T1D", trait="Type 1 diabetes", locus=paste0("T1D::", `Signal name`),
  chr=paste0("chr", as.character(`Chrom.`)), pos=as.integer(`Position (hg38)`),
  pip=as.numeric(PPA), consequence=NA_character_, variant_id=as.character(`Unique ID (hg38)`)
)]
t1_cluster <- cluster_signal_loci(t1)
t1_obj <- apply_qc(t1_cluster$variants, t1_cluster$loci, FALSE)
write_reference("T1D", t1_obj,
  "Chiou et al. 2021 Supplementary Data 1, 99% credible sets for 136 signals",
  "GRCh38 supplied", "+/-500 kb around highest-PIP variant per signal")

message("Preparing inflammatory bowel disease")
ib_file <- file.path(ref_root,"IBD","Huang_2017_supplementary_table1.xlsx")
ib0 <- as.data.table(read_excel(ib_file, sheet="list of variants"))
ib <- ib0[, .(
  source="IBD", trait="Inflammatory bowel disease", locus=paste0("IBD::",HD),
  chr=paste0("chr",chr), pos=as.integer(position), pip=as.numeric(P_mean_95),
  consequence=paste(Func.refGene, ExonicFunc.refGene, sep=";"), variant_id=as.character(variant)
)]
ib <- ib[!is.na(pos) & !is.na(pip)]
ib_lift <- lift_points(ib)
ib_lift_n_input <- attr(ib_lift,"lift_n_input")
ib_lift_n_output <- attr(ib_lift,"lift_n_output")
ib_cs <- as.data.table(read_excel(ib_file, sheet="list of credible sets"))
ib_reg <- tstrsplit(as.character(ib_cs$`region (mb)`), "-", fixed=TRUE)
ib_loci0 <- ib_cs[, .(source="IBD", trait="Inflammatory bowel disease",
                      locus=paste0("IBD::",HD), chr=paste0("chr",chr),
                      start=as.integer(round(as.numeric(ib_reg[[1]])*1e6)),
                      end=as.integer(round(as.numeric(ib_reg[[2]])*1e6)))]
ib_loci0 <- unique(ib_loci0, by="locus")
ib_loci0 <- ib_loci0[!is.na(start) & !is.na(end) & start <= end &
                     grepl("^chr([1-9]|1[0-9]|2[0-2]|X)$", chr)]
gr <- GRanges(ib_loci0$chr, IRanges(ib_loci0$start, ib_loci0$end))
lo <- liftOver(gr, chain19to38)
ib_loci <- rbindlist(lapply(seq_along(lo), function(i) {
  z <- lo[[i]]
  if (!length(z)) return(NULL)
  data.table(source=ib_loci0$source[i], trait=ib_loci0$trait[i],
             locus=ib_loci0$locus[i], chr=as.character(seqnames(z))[1],
             start=min(start(z)), end=max(end(z)))
}))
ib_obj <- apply_qc(ib_lift, ib_loci, TRUE)
write_reference("IBD", ib_obj,
  "Huang et al. 2017 Supplementary Table 1 variant and credible-set tables",
  "GRCh37/hg19", "published fine-mapping region_start/region_end",
  ib_lift_n_input, ib_lift_n_output)

all_manifest <- rbindlist(lapply(c("FinnGen","RA","T1D","IBD"), function(s)
  fread(file.path(ref_root,s,"processed","reference_manifest.tsv"))), fill=TRUE)
fwrite(all_manifest, file.path(ref_root,"SCENT_style_reference_manifest.tsv"), sep="\t")
message("Wrote standardised references and manifests under ", ref_root)
