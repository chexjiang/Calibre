#!/usr/bin/env Rscript

suppressPackageStartupMessages({library(data.table);library(ggplot2)})
args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(if(length(args)) args[[1]] else ".",mustWork=TRUE)
source <- if(length(args)>=2) args[[2]] else "RA"
base <- file.path(root,"results","validation","GWAS","SCENT_style",source)
figdir <- file.path(base,"figures");dir.create(figdir,recursive=TRUE,showWarnings=FALSE)
x <- fread(file.path(base,"tables","scent_summary.tsv"))
cols <- c("EpiLink-Poisson"="#E47159","EpiLink-NB"="#C87878",
          "scMultiMap"="#edd089","Signac (corrected)"="#54b2a2",
          "Signac (uncorrected)"="#78A878")
lev <- names(cols)
x[,method_label:=factor(method_label,levels=lev)]
x[,dataset:=factor(dataset,levels=unique(dataset))]

theme_bench <- theme_bw(base_size=10)+theme(
  panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),
  axis.text.x=element_text(angle=35,hjust=1,vjust=1),
  strip.background=element_rect(fill="#F2F2F2",colour="#777777"),
  legend.position="none",plot.title=element_text(face="bold"))

draw <- function(z,title) ggplot(z,aes(method_label,display_score,color=method_label))+
  geom_hline(yintercept=1,linetype=3,color="grey70",linewidth=.35)+
  geom_errorbar(aes(ymin=display_ci_low,ymax=display_ci_high),width=.16,linewidth=.45)+
  geom_point(size=2.3)+scale_color_manual(values=cols,drop=FALSE)+
  labs(x=NULL,y=expression(log[2](1+enrichment)),title=title)+theme_bench

p <- draw(x[scope=="pooled"],paste0(source," SCENT-style enrichment - pooled cell types"))+
  facet_wrap(~dataset,scales="free_y",ncol=3)
ggsave(file.path(figdir,"pooled_method_comparison.pdf"),p,width=11,height=max(7,3.0*ceiling(uniqueN(x[scope=="pooled"]$dataset)/3)))

pdf(file.path(figdir,"celltype_method_comparison.pdf"),width=11,height=8.5,onefile=TRUE)
for(d in unique(as.character(x[scope=="celltype"]$dataset))) {
  z <- x[scope=="celltype" & as.character(dataset)==d]
  p <- draw(z,paste0(source," SCENT-style enrichment - ",d))+facet_wrap(~celltype,scales="free_y",ncol=3)
  print(p)
}
dev.off()

# A compact cross-dataset overview without confidence bars. This is secondary
# to the count table and point/error-bar figures, but makes broad patterns easy
# to scan without forcing a shared raw-enrichment axis.
hp <- ggplot(x[scope=="pooled"],aes(method_label,dataset,fill=display_score))+
  geom_tile(color="white",linewidth=.5)+
  geom_text(aes(label=sprintf("%.2f",display_score)),size=2.7)+
  scale_fill_gradient2(low="#F7FBFF",mid="#FFF7BC",high="#B2182B",midpoint=1)+
  labs(x=NULL,y=NULL,fill=expression(log[2](1+E)),title=paste0(source," pooled enrichment overview"))+
  theme_minimal(base_size=10)+theme(axis.text.x=element_text(angle=35,hjust=1),
    panel.grid=element_blank(),plot.title=element_text(face="bold"))
ggsave(file.path(figdir,"pooled_heatmap.pdf"),hp,width=8.5,height=max(5,0.38*uniqueN(x[scope=="pooled"]$dataset)+2))

ranked <- x[,.(rank=frank(-overall_enrichment,ties.method="average"),overall_enrichment,
                ci_low,ci_high,display_score),by=.(source,scope,dataset,celltype,method_label)]
fwrite(ranked,file.path(base,"tables","method_ranks.tsv"),sep="\t")
message("Wrote PDF figures to ",figdir)
