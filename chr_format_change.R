# function for chromosome format change

N2chr<-function(data){
    data[,chr:=as.character(chr)]
    data[chr=="1"]$chr="chr1"
    data[chr=="2"]$chr="chr2"
    data[chr=="3"]$chr="chr3"
    data[chr=="4"]$chr="chr4"
    data[chr=="5"]$chr="chr5"
    data[chr=="Mt"]$chr="chrM"
    data[chr=="Pt"]$chr="chrC"
    return(data)
}
chr2N<-function(data){
    data[chr=="chr1"]$chr="1"
    data[chr=="chr2"]$chr="2"
    data[chr=="chr3"]$chr="3"
    data[chr=="chr4"]$chr="4"
    data[chr=="chr5"]$chr="5"
    data[chr=="chrM"]$chr="Mt"
    data[chr=="chrC"]$chr="Pt"
    return(data)
}

chr2Chr<-function(data){
    data[,chr:=str_c("C",str_sub(chr,2,4))]
    return(data)
}

Chr2chr<-function(data){
    data[,chr:=str_c("c",str_sub(chr,2,4))]
    return(data)
}

N2Chr<-function(data){
    return(chr2Chr(N2chr(data)))
}
Chr2N<-function(data){
    return(chr2N(Chr2chr(data)))
}

