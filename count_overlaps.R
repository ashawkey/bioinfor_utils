require(GenomicRanges)
#` @gquery: chip-peaks
#' @gsubject: gff 
findOverlaps_collapsed<-function(gquery, gsubject){
    total <- length(gquery)
    as.data.table(findOverlaps(gquery,gsubject)) -> ops
    gquery[ops$queryHits] -> gq_expand
    gq_expand$overlap <- gsubject[ops$subjectHits]$type
    gq_expand <- as.data.table(gq_expand)
    melt(gq_expand,measure.vars="overlap") -> gq_expand
    dcast(gq_expand,seqnames+start+end+width+strand~variable,fun=paste,collapse="|") -> gq_ann
    # make a table of statistics
    type=c("gene","transposable_element","transposable_element_gene","transposable_element_within_gene","pseudogene")   
    res=data.table(type=type,percent=0)
    for(t in type){
        nrow(gq_ann[str_detect(overlap,t)])/total -> res[type==t,"percent"]
    }
    res=rbind(res,data.table(type="intergenic", percent=(total-nrow(gq_ann))/total))
    res$type=c("gene","TE out gene","TEG","TE in gene","pseudogene","intergenic")
    res$type=factor(res$type,levels=c("gene","TEG","TE out gene","TE in gene","pseudogene","intergenic"))
    return(res)
}
