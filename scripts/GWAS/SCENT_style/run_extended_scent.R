#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
})

args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(if (length(args)) args[[1]] else ".",mustWork=TRUE)
source <- if (length(args)>=2) args[[2]] else "RA"
n_boot <- if (length(args)>=3) as.integer(args[[3]]) else 1000L
seed <- if (length(args)>=4) as.integer(args[[4]]) else 20261003L
stopifnot(source %chin% c("FinnGen","RA","T1D","IBD"))
set.seed(seed)

gwas <- file.path(root,"results","validation","GWAS")
base <- file.path(gwas,"SCENT_style",source)
dir.create(file.path(base,"results"),recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(base,"tables"),recursive=TRUE,showWarnings=FALSE)

manifest <- fread(file.path(gwas,"SCENT_style","config","extended_annotation_manifest.tsv"))
ref <- readRDS(file.path(gwas,"reference",source,"processed","scent_reference.rds"))
causal <- unique(ref$causal[,.(trait,locus,chr,pos)],by=c("trait","locus","chr","pos"))
loci <- unique(ref$loci[,.(trait,locus,chr,start,end)])
common <- readRDS(file.path(gwas,"reference","1000G","1000G_EUR_MAFgt0.01.hg38.rds"))
common <- unique(common[,.(chr=as.character(chr),pos=as.integer(pos))],by=c("chr","pos"))
common_gr <- GRanges(common$chr,IRanges(common$pos,common$pos))

read_bed <- function(path,id) {
  x <- fread(path,header=FALSE,col.names=c("chr","start0","end"))
  if (!nrow(x)) return(data.table(chr=character(),start=integer(),end=integer(),annotation_id=character()))
  x[,.(chr,start=as.integer(start0)+1L,end=as.integer(end),annotation_id=id)]
}
ann_dt <- rbindlist(lapply(seq_len(nrow(manifest)),function(i)
  read_bed(manifest$bed_path[i],manifest$annotation_id[i])))
ann_gr <- GRanges(ann_dt$chr,IRanges(ann_dt$start,ann_dt$end),annotation_id=ann_dt$annotation_id)
ann_ids <- manifest$annotation_id

count_annotations <- function(points) {
  ans <- setNames(integer(length(ann_ids)),ann_ids)
  if (!nrow(points)) return(ans)
  pgr <- GRanges(points$chr,IRanges(points$pos,points$pos))
  h <- findOverlaps(pgr,ann_gr,ignore.strand=TRUE)
  if (length(h)) {
    # BED intervals were reduced within each annotation; nevertheless count a
    # variant only once per annotation to make the invariant explicit.
    pairs <- unique(data.table(point=queryHits(h),id=mcols(ann_gr)$annotation_id[subjectHits(h)]))
    tt <- table(pairs$id); ans[names(tt)] <- as.integer(tt)
  }
  ans
}

common_in_intervals <- function(intervals) {
  if (!nrow(intervals)) return(common[0])
  igr <- reduce(GRanges(intervals$chr,IRanges(intervals$start,intervals$end)),ignore.strand=TRUE)
  h <- findOverlaps(common_gr,igr,ignore.strand=TRUE)
  common[unique(queryHits(h))]
}

traits <- sort(intersect(unique(causal$trait),unique(loci$trait)))
rows <- vector("list",length(traits))
for (ti in seq_along(traits)) {
  tr <- traits[ti]
  ca <- unique(causal[trait==tr,.(chr,pos)])
  bg <- common_in_intervals(loci[trait==tr])
  nca <- nrow(ca); nbg <- nrow(bg)
  ca_n <- count_annotations(ca); bg_n <- count_annotations(bg)
  z <- data.table(annotation_id=ann_ids,n_causal_total=nca,n_causal_in=unname(ca_n),
                  n_common_total=nbg,n_common_in=unname(bg_n))
  z[, enrichment := fifelse(n_causal_total>0 & n_common_total>0 & n_common_in>0,
    (n_causal_in/n_common_in)/(n_causal_total/n_common_total),NA_real_)]
  z[, `:=`(source=source,trait=tr)]
  rows[[ti]] <- z
  if (ti %% 25L == 0L || ti==length(traits))
    message(source,": trait ",ti,"/",length(traits),"; causal=",nca,"; common=",nbg)
}
per_trait <- rbindlist(rows)
per_trait <- merge(per_trait,manifest,by="annotation_id",all.x=TRUE)
fwrite(per_trait,file.path(base,"results","per_trait_enrichment.tsv.gz"),sep="\t")

if (length(traits)>1L) {
  summary <- per_trait[is.finite(enrichment), {
    vals <- enrichment
    b <- replicate(n_boot,mean(sample(vals,length(vals),replace=TRUE)))
    list(n_bootstrap_units=length(vals),overall_enrichment=mean(vals),
         median_trait_enrichment=median(vals),
         ci_low=unname(quantile(b,.025)),ci_high=unname(quantile(b,.975)),
         n_causal_in=sum(n_causal_in),n_causal_total=sum(n_causal_total),
         n_common_in=sum(n_common_in),n_common_total=sum(n_common_total))
  },by=.(annotation_id,dataset,celltype,scope,method,method_label,correction)]
  summary[, bootstrap_unit:="trait"]
  per_locus <- data.table()
} else {
  # One disease is one trait. Preserve the SCENT pooled count ratio, and use
  # loci rather than traits as the resampling unit for uncertainty.
  lr <- vector("list",nrow(loci))
  for (li in seq_len(nrow(loci))) {
    loc <- loci[li]
    ca <- unique(causal[locus==loc$locus,.(chr,pos)])
    bg <- common_in_intervals(loc)
    ca_n <- count_annotations(ca); bg_n <- count_annotations(bg)
    lr[[li]] <- data.table(locus=loc$locus,annotation_id=ann_ids,
      n_causal_total=nrow(ca),n_causal_in=unname(ca_n),
      n_common_total=nrow(bg),n_common_in=unname(bg_n))
    if (li %% 20L==0L || li==nrow(loci)) message(source,": locus ",li,"/",nrow(loci))
  }
  per_locus <- rbindlist(lr)
  fwrite(per_locus,file.path(base,"results","per_locus_counts.tsv.gz"),sep="\t")
  calc_ratio <- function(z) {
    if (sum(z$n_common_in)==0 || sum(z$n_causal_total)==0) return(NA_real_)
    (sum(z$n_causal_in)/sum(z$n_common_in))/(sum(z$n_causal_total)/sum(z$n_common_total))
  }
  summary <- per_locus[, {
    z <- .SD
    bs <- replicate(n_boot,{
      pick <- sample(seq_len(.N),.N,replace=TRUE)
      calc_ratio(z[pick])
    })
    list(n_bootstrap_units=.N,overall_enrichment=calc_ratio(z),
         median_trait_enrichment=NA_real_,ci_low=unname(quantile(bs,.025,na.rm=TRUE)),
         ci_high=unname(quantile(bs,.975,na.rm=TRUE)),
         n_causal_in=sum(n_causal_in),n_causal_total=sum(n_causal_total),
         n_common_in=sum(n_common_in),n_common_total=sum(n_common_total))
  },by=annotation_id]
  summary <- merge(summary,manifest[,.(annotation_id,dataset,celltype,scope,method,
                                        method_label,correction)],by="annotation_id")
  summary[, bootstrap_unit:="locus"]
}
summary[, `:=`(source=source,n_bootstrap=n_boot,seed=seed,
              display_score=log2(1+overall_enrichment),
              display_ci_low=log2(1+ci_low),display_ci_high=log2(1+ci_high))]
setcolorder(summary,c("source","scope","dataset","celltype","method","method_label",
                      "correction","overall_enrichment","ci_low","ci_high","display_score",
                      "display_ci_low","display_ci_high"))
fwrite(summary,file.path(base,"tables","scent_summary.tsv"),sep="\t")
fwrite(summary[,.(source,scope,dataset,celltype,method_label,correction,
                  n_causal_in,n_causal_total,n_common_in,n_common_total,
                  overall_enrichment,ci_low,ci_high)],
       file.path(base,"tables","count_and_enrichment_table.tsv"),sep="\t")
fwrite(data.table(source=source,n_traits=length(traits),n_loci=nrow(loci),
                  n_causal=nrow(causal),n_annotations=nrow(manifest),
                  bootstrap_unit=unique(summary$bootstrap_unit),n_bootstrap=n_boot),
       file.path(base,"results","run_audit.tsv"),sep="\t")
message("Finished ",source,": ",nrow(summary)," annotation summaries")
