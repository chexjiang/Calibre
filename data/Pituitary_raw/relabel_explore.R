suppressPackageStartupMessages({library(Matrix);library(ggplot2);library(irlba);library(uwot)})
out_dir <- file.path("data","Pituitary")
rna <- readRDS(file.path(out_dir,"Pituitary.rna.rds")); meta <- readRDS(file.path(out_dir,"Pituitary.meta.rds"))
keep <- intersect(rownames(meta),colnames(rna)); rna <- rna[,keep,drop=FALSE]; meta <- meta[keep,,drop=FALSE]
lib <- Matrix::colSums(rna); norm <- rna %*% Diagonal(x=1e4/pmax(lib,1)); norm@x <- log1p(norm@x)
sets <- list(Somatotropes=c("GH1","GH2","GHRHR"),Lactotropes=c("PRL","PIT1","POU1F1"),Corticotropes=c("POMC","TBX19"),Gonadotropes=c("LHB","FSHB","GNRHR"),Thyrotropes=c("TSHB","TRHR"))
scores <- sapply(sets,function(gs){g<-intersect(gs,rownames(norm));if(!length(g))return(rep(0,ncol(norm)));Matrix::colMeans(norm[g,,drop=FALSE])})
meta$newCT <- colnames(scores)[max.col(scores,ties.method="first")]; meta$celltype <- factor(meta$newCT,levels=c("Somatotropes","Corticotropes","Lactotropes","Thyrotropes","Gonadotropes")); meta$newCT <- as.character(meta$celltype)
keep_types <- names(table(meta$celltype))[table(meta$celltype)>500]; keep <- rownames(meta)[meta$celltype %in% keep_types]; meta <- meta[keep,,drop=FALSE]; rna <- rna[,keep,drop=FALSE]
saveRDS(meta,file.path(out_dir,"Pituitary.meta.rds"),compress="gzip")
ct <- as.data.frame(table(meta$newCT),stringsAsFactors=FALSE); colnames(ct)<-c("celltype","n_cells"); write.table(ct,file.path(out_dir,"Pituitary.celltype_counts.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
lib <- Matrix::colSums(rna); z <- rna %*% Diagonal(x=1e4/pmax(lib,1)); z@x<-log1p(z@x); va<-Matrix::rowMeans(z^2)-Matrix::rowMeans(z)^2; hv<-order(va,decreasing=TRUE)[seq_len(min(2000,length(va)))]
pc<-irlba::prcomp_irlba(t(z[hv,,drop=FALSE]),n=20,center=TRUE,scale.=FALSE); um<-uwot::umap(pc$x,n_neighbors=30,min_dist=.3,metric="cosine",n_threads=8,verbose=FALSE); emb<-as.data.frame(um); colnames(emb)<-c("UMAP1","UMAP2"); emb$celltype<-meta$celltype[match(rownames(emb),rownames(meta))]; nb<-table(emb$celltype); emb$label<-paste0(emb$celltype," (",nb[emb$celltype],")")
p<-ggplot(emb,aes(UMAP1,UMAP2,colour=label))+geom_point(size=.28,alpha=.5)+theme_classic(base_size=11)+labs(title="Pituitary",subtitle="RNA PCA/UMAP after count QC",colour=NULL)+guides(colour=guide_legend(override.aes=list(size=3,alpha=1))); ggsave(file.path(out_dir,"Pituitary.celltype_umap.png"),p,width=6.5,height=5,dpi=220)
cat("Relabeled pituitary:",nrow(meta),"cells\n"); print(ct)
