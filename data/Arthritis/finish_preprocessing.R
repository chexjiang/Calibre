# Finish the arthritis outputs after the large H5AD-to-RDS conversion.
suppressPackageStartupMessages({
  library(Matrix); library(Seurat); library(Signac); library(ggplot2); library(rhdf5)
})
source(file.path("R", "helper.R"))
out_dir <- file.path("data", "Arthritis")
rna <- readRDS(file.path(out_dir, "Arthritis.rna.rds"))
meta <- readRDS(file.path(out_dir, "Arthritis.meta.rds"))
annotation <- readRDS(file.path("data", "hg38_EnsDb_Hsapiens_v86_annotation.rds"))
obj <- SeuratObject::CreateSeuratObject(counts=rna, assay="RNA", meta.data=meta)
peak_names <- as.character(rhdf5::h5read(
  file.path(out_dir, "GSE243917_peaksXcells_ALL_raw.h5ad"), "var/_index"))
peaks <- Signac::StringToGRanges(peak_names, sep=c(":","-"))
gene_range <- CollapseToLongestTranscript(annotation)
gene_window <- Signac::Extend(gene_range, upstream=5e5, downstream=5e5)
common_chr <- intersect(GenomeInfoDb::seqlevels(peaks), GenomeInfoDb::seqlevels(gene_window))
peaks <- GenomeInfoDb::keepSeqlevels(peaks, common_chr, pruning.mode="coarse")
gene_window <- GenomeInfoDb::keepSeqlevels(gene_window, common_chr, pruning.mode="coarse")
hits <- GenomicRanges::findOverlaps(peaks, gene_window)
gene_peak_df <- data.frame(
  gene=gene_window$gene_name[S4Vectors::subjectHits(hits)],
  peak=peak_names[S4Vectors::queryHits(hits)], stringsAsFactors=FALSE)
gene_peak_df <- gene_peak_df[gene_peak_df$gene %in% rownames(rna),,drop=FALSE]
rownames(gene_peak_df) <- NULL
gene_peak_df <- gene_peak_df[!duplicated(gene_peak_df),c("gene","peak"),drop=FALSE]
gene_peak_df$peak <- CanonicalPeakNames(gene_peak_df$peak)
saveRDS(gene_peak_df,file.path(out_dir,"Arthritis.gene_peak_df.rds"),compress="gzip")
desired_order <- c("T","Myeloid","B","Fibroblast","Endothelial")
counts <- as.data.frame(table(meta$newCT),stringsAsFactors=FALSE)
colnames(counts) <- c("celltype","n_cells")
counts <- counts[order(match(counts$celltype,desired_order)),,drop=FALSE]
write.table(counts,file.path(out_dir,"Arthritis.celltype_counts.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
qc_metrics <- c("nCount_RNA","nFeature_RNA","nCount_ATAC","nFeature_ATAC")
qc_summary <- do.call(rbind,lapply(qc_metrics,function(nm){x<-as.numeric(meta[[nm]]);data.frame(metric=nm,n=length(x),min=min(x),median=median(x),mean=mean(x),max=max(x))}))
write.table(qc_summary,file.path(out_dir,"Arthritis.qc_summary.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
qc_thresholds <- data.frame(dataset="arthritis-tissue",tissue="joint",RNA_count="> 500 genes",ATAC_count="> 5000",nucleosome="10,000 reads (paper table)",TSS_enrichment=">= 10% of reads in promoters",source_matrix="GEO post-QC matrices; no second QC pass")
write.table(qc_thresholds,file.path(out_dir,"Arthritis.qc_thresholds.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
explore <- NormalizeData(obj,verbose=FALSE)
explore <- FindVariableFeatures(explore,nfeatures=2000,verbose=FALSE)
explore <- ScaleData(explore,features=VariableFeatures(explore),verbose=FALSE)
explore <- RunPCA(explore,features=VariableFeatures(explore),npcs=30,verbose=FALSE)
explore <- RunUMAP(explore,dims=1:20,n.neighbors=30,min.dist=.3,verbose=FALSE)
emb <- as.data.frame(Embeddings(explore,"umap")[,1:2,drop=FALSE])
colnames(emb) <- c("UMAP1","UMAP2")
emb$celltype <- meta$celltype[match(rownames(emb),rownames(meta))]
n_by_type <- table(emb$celltype)
emb$label <- paste0(emb$celltype," (",n_by_type[emb$celltype],")")
p <- ggplot(emb,aes(UMAP1,UMAP2,colour=label))+geom_point(size=.22,alpha=.45)+
  theme_classic(base_size=11)+labs(title="Arthritis tissue (joint)",
  subtitle="RNA PCA/UMAP; GEO post-QC shared cells",colour=NULL)+
  guides(colour=guide_legend(override.aes=list(size=3,alpha=1)))
ggsave(file.path(out_dir,"Arthritis.celltype_umap.png"),p,width=6.5,height=5,dpi=220)
qc_long <- do.call(rbind,lapply(qc_metrics,function(nm)data.frame(celltype=meta$celltype,metric=nm,value=as.numeric(meta[[nm]]))) )
q <- ggplot(qc_long,aes(celltype,value,fill=celltype))+geom_violin(scale="width",trim=TRUE)+
  facet_wrap(~metric,scales="free_y",ncol=2)+theme_classic(base_size=10)+
  theme(axis.text.x=element_text(angle=35,hjust=1),legend.position="none")+
  labs(title="Arthritis data exploration: count distributions",x=NULL,y=NULL)
ggsave(file.path(out_dir,"Arthritis.qc_violin.png"),q,width=9,height=6,dpi=220)
saveRDS(explore,file.path(out_dir,"Arthritis.exploration_seurat.rds"),compress="gzip")
writeLines(c("Arthritis tissue benchmark preparation","",
  "Source: GEO GSE243917. Complete raw gene/peak H5AD matrices are already post-QC.",
  "Thresholds are in Arthritis.qc_thresholds.tsv.",
  "The complete GEO object has 654 mural cells; these are excluded to reproduce the published 30,893-cell SCENT panel. Bcell+plasma and Tcell+NK are merged.",
  "ATAC.rds is a peak-by-cell count matrix; SCENT should binarize it method-side.",
  "Candidate gene_peak_df uses hg38 gene bodies +/-500 kb and the shared ATAC peak set."),
  file.path(out_dir,"README.md"))
cat("Finished arthritis:",ncol(rna),"cells,",length(peak_names),"peaks,",nrow(gene_peak_df),"candidate pairs.\n")
