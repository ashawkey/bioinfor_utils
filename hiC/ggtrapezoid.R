# some utils for coord change ---

# chrpos(string obj): chr:st,ed
# extension: chr:st,chr:ed
id2chrpos<-function(id, bed){
    colnames(bed)=c("chr","start","end","idx")
    return(str_c(bed[idx==id,chr],":",bed[idx==id,start],",",bed[idx==id,end]))
}

region2id<-function(region, bed){
    colnames(bed)=c("chr","start","end","idx")
    interval=bed[1,end]-bed[1,start]
    Chr=str_replace(region,"(.*):(.*),(.*)","\\1")
    Start=str_replace(region,"(.*):(.*),(.*)","\\2")
    End=str_replace(region,"(.*):(.*),(.*)","\\3")
    
    if(Start==Chr & End==Chr){
        Chr=str_sub(Chr,1,4)
        Start=bed[chr==Chr]$start[1]
        End=bed[chr==Chr]$end[length(bed[chr==Chr]$end)]
        id_start=bed[chr==Chr]$idx[1]
        id_end=bed[chr==Chr]$idx[length(bed[chr==Chr]$idx)]
    } else {
        Start=as.integer(Start)
        End=as.integer(End)
        tmp = bed[chr==Chr&start<=End&end>=Start]$idx
        id_start=tmp[1]
        id_end=tmp[length(tmp)]
    }
    return(c(id_start,id_end,Chr,Start,End))
}

# region is a chrpos string
ggtrapezoid2<-function(mat,
                       bed,
                       region="all",
                       dist=-1,
                       type="full",
                       color="red",
                       by=50,
                       transposed=F,
                       collision_mininal_distance=10
                       ){
    # input transform
    colnames(bed) = c("chrx","startx","endx","idx")
    if(region=="all"){
        cat("[INFO] Region is set to whole genome as default.\n")
        Chr="Genome"
        start = bed$idx[1]
        end = bed$idx[length(bed$idx)]
    }else{
        tmp = region2id(region,bed)
        start = as.integer(tmp[1])
        end = as.integer(tmp[2])
        Chr = tmp[3]
        Start = as.integer(tmp[4])
        End = as.integer(tmp[5])
    }
    # constants
    binSize=bed[1,endx]-bed[1,startx]
    by=by*binSize
    diag_len = end - start + 1
    total_num = diag_len*diag_len
    if(dist == -1 | dist > diag_len) {
        dist=diag_len
        cat("[INFO] Dist is set to max value ",diag_len," as default.\n")
    }
    # error input 
    if(diag_len > nrow(bed))
    {
        return("[ERROR] Region is larger than input data range.\n")
    }
    # calculate all centers of diamond polygons to be plotted
    # use vectorized operation to speed up
    id=1:total_num
    i=ceiling(id/diag_len)
    j=id-diag_len*(i-1)
    # i is col, j is row, rotate 45 degree
    diamonds_center=data.table(id=id,x=i+j-1,y=i-j)
    
    # from center to calculate the vertices
    diamonds_up=diamonds_center
    diamonds_up$y=diamonds_up$y+1
    diamonds_right=diamonds_center
    diamonds_right$x=diamonds_right$x+1
    diamonds_down=diamonds_center
    diamonds_down$y=diamonds_down$y-1
    diamonds_left=diamonds_center
    diamonds_left$x=diamonds_left$x-1

    diamonds_vertices=rbind(diamonds_up,diamonds_right,diamonds_down,diamonds_left)

    # crop the region interested
    colnames(mat)=c("i","j","val")
    mmat=mat[i<=end&i>=start&j<=end&j>=start]
    mmat[,i:=i-start+1]
    mmat[,j:=j-start+1]
    if(!transposed){
        mmat[,id:=(j-1)*diag_len+i]
    } else mmat[,id:=(i-1)*diag_len+j]
    
    mmat[diamonds_vertices,on="id"][,c(4:6,3)]->datapoly
    datapoly[is.na(val)]$val=0

    # log transformation
    datapoly[,val:=log10(val+1)]
    
    # merge to create plot data
    datapoly=datapoly[abs(y)<=dist]
    if(type=="up") datapoly=datapoly[y>=0]
    else if(type=="bottum") datapoly=datapoly[y<=0]
    
    # ggplot2
    xmax=2*diag_len-1
    tmp=bed[startx%%by==0 & idx>=start & idx<=end]
    # transformed max id of current plot
    idmax=diag_len-1+tmp$idx[1]
    tmp=rbind(tmp,bed[idx==idmax])
    # modification to avoid ugly overlapping labels
    tmp$nxt=tmp$idx[c(2:length(tmp$idx),1)]
    tmp[,dis:=nxt-idx]
    tmp=tmp[dis<0 | dis>collision_mininal_distance]
    # transformation from idx to real x coordinates.
    myBreaks=2*(tmp$idx-tmp$idx[1])
    myLabels=round(tmp$start/(binSize))

    p=ggplot(datapoly,aes(x,y))+
        geom_polygon(aes(fill=val,group=id))+
        scale_x_continuous(position="top",
                           limits=c(0,xmax),
                           expand=c(0,0),
                           breaks=myBreaks,
                           labels=myLabels
                           )+
        
        xlab(Chr)+
        ylim(min(datapoly$y),max(datapoly$y))+
        scale_fill_gradientn(colours=c("white",color))+
        coord_fixed(ratio=1)+
        theme(panel.grid=element_blank(),
              panel.background=element_blank(),
              axis.text.y=element_blank(),
              axis.title.y=element_blank(),
              axis.ticks.y=element_blank(),
              axis.line.x = element_line(size=0.7,lineend="round"),
              axis.ticks.length=unit(.25,'cm')
              )
    return(p)
}
