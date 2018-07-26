#! /usr/bin/env Rscript

Args<-commandArgs(trailingOnly = T)

if (length(Args)!=3) {
    stop("Usage:Rscript enrich.R totalRegion.rds interestedRegion.rds binData.bedGraph outputName", call.=F)
}

## Import necessary packages ----
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(Rcpp))

sourceCpp("sumBinOverRegion.cc")

# -----------------------
# the input's chrM/Mt/ChrM must be the same.
# interestedRegions should provide "cls" column for plotting
# ----------------------

readRDS(Args[1])->regions
readRDS(Args[2])->interestedRegions
fread(Args[3])->binData

colnames(binData)=c("seqnames","start","end","val")
as.data.table(regions)->regions
as.data.table(interestedRegions)->interestedRegions

regions$cls="genome"
regions$val=sumBinOverRegion(regions,binData)
regions[,val:=val/width]

interestedRegions$val=sumBinOverRegion(interestedRegions,binData)
interestedRegions[,val:=val/width]

rbind(regions,interestedRegions)[,c("cls","val")]->pdata

ggplot(pdata,aes(as.factor(cls),val,fill=as.factor(cls)))+
    geom_boxplot()+
    geom_jitter(data=pdata[cls!="genome"],position=position_jitter(0.2),alpha=0.5,size=3)+
    scale_y_log10()+
    xlab("class")+
    ylab("density")+
    ggtitle(str_c(outputName," Enrichment"))+
    theme_bw()+
    theme(legend.position="none",text=element_text(size=16),plot.title=element_text(hjust=0.5))+
    geom_signif(comparisons=list(c("genome","1"),c("genome","2"),c("genome","3"),c("genome","4"),c("genome","5")),map_signif_level=T,step_increase=0.1)

ggsave(str_c(outputName,".png"))


