# transform any bedGraph file to uniformly binned bedGraph ---

genome_tiler<-function(chrsize,bin){
    chrsize[,V3:=1]
    chrsize[,c(1,3,2)]->chrsize
    colnames(chrsize)=c("chr","start","end")
    makeGRangesFromDataFrame(chrsize)->gchr
    unlist(tile(gchr,width=bin))->gchr_tiles
    return(gchr_tiles)
}

bedGraph_transform<-function(bg1,bg2,fun="mean"){
    as.data.table(findOverlapPairs(bg2,bg1))->tmp
    tmp[,c(1,2,3,5,7,8,9,12)]->tmp
    colnames(tmp)=c("chr","start","end","strand","st2","ed2","wid2","val")
    melt(tmp,measure="val")->tmp
    if(fun=="mean"){
        fun=eval(parse(text=fun))
        tmp[,lw:=as.integer(ed2-start+1)]
        tmp[,rw:=as.integer(end-st2+1)]
        tmp[wid2>lw,wid2:=lw]
        tmp[wid2>rw,wid2:=rw]
        tmp[,value:=value*wid2]
        dcast(tmp,chr+start+end+strand~variable,fun=sum)->res
        res[,width:=end-start+1]
        res[,val:=val/width]
    } else {
        fun=eval(parse(text=fun))
        dcast(tmp,chr+start+end+strand~variable,fun=fun)->res
    }
    return(res)
}

getLog2<-function(res1,res2){
    res1[,id:=str_c(chr,":",start,"-",end)]
    res2[,id:=str_c(chr,":",start,"-",end)]
    res1[res2,on="id"]->res
    res[,c(1,2,3,5,12)]->res
    res[,lgTI:=log2(val/i.val)]
    return(res)
}

enrich_preprocess<-function(lst, bedgraph){
    colnames(bedgraph)=c("chr","start","end","val")
    pdata=data.table()
    for(i in 1:length(lst)){ 
        colnames(lst[[i]])=c("chr","start","end")
        makeGRangesFromDataFrame(lst[[i]])->bg2
        makeGRangesFromDataFrame(bedgraph,keep.extra.columns=T)->bg1
        bedGraph_transform(bg1,bg2)->res
        res$type=names(lst)[i]
        pdata=rbind(pdata,res)
    }
    return(pdata)
}

