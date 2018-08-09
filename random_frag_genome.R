random_frag<-function(chrSize, num, len, seed=1){
    set.seed(seed)
    chrSize <- as.data.table(chrSize)
    colnames(chrSize)=c("chr","size")
    totalBases=sum(chrSize$size)
    chrSize[,fragNum:=as.integer((size/totalBases)*num)]
    bed=data.table()
    for(i in 1:nrow(chrSize)){
        if(chrSize[i,fragNum]==0) next
        tmp=data.table(chr=rep(chrSize[i,chr],chrSize[i,fragNum]))
        tmp[,start:=as.integer(runif(chrSize[i,fragNum], 0, chrSize[i,size]-len))]
        tmp[,end:=start+len]
        tmp[,strand:="*"]
        bed=rbind(bed,tmp)
    }
    return(bed[order(chr,start,end)])
}
