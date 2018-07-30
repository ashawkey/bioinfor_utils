# some utils ---

id2chrpos<-function(id, bed){
    colnames(bed)=c("chr","start","end","idx")
    return(str_c(bed[idx==id,chr],":",bed[idx==id,start],",",bed[idx==id,end]))
}

#' region defination:
#' :                 whole genome
#' chrX-chrY         whole chrX to chrY
#' chrX              whole chrX
#' chrX:aaaaa-bbbbb  chrX subregion aaaaa to bbbbb
region2id<-function(region, bed){
    colnames(bed)=c("chr","start","end","idx")
    interval=bed[1,end]-bed[1,start]
    if(region==":"){
        Chr="genome"
        id_start = bed$idx[1]
        id_end = bed$idx[length(bed$idx)]
    } else if(str_detect(region,"-")){
        if(str_detect(region,":")){
            Chr=str_replace(region,"(.*):(.*)-(.*)","\\1")
            Start=as.integer(str_replace(region,"(.*):(.*)-(.*)","\\2"))
            End=as.integer(str_replace(region,"(.*):(.*)-(.*)","\\3"))
            tmp = bed[chr==Chr&start<=End&end>=Start]$idx
            id_start=tmp[1]
            id_end=tmp[length(tmp)]
        }else{
            Chr=region
            Chr1=str_replace(region,"(.*)-(.*)","\\1")
            Chr2=str_replace(region,"(.*)-(.*)","\\2")
            id_start=bed[chr==Chr1]$idx[1]
            id_end=bed[chr==Chr2]$idx[length(bed[chr==Chr2]$idx)]
        }
    } else {
        Chr=region
        id_start = bed[chr==Chr]$idx[1]
        id_end = bed[chr==Chr]$idx[length(bed[chr==Chr]$idx)]
    }
    return(c(id_start,id_end,Chr))
}

#' ggtrapezoid
#' @mat: sparse matrix format
#' @bed: bedGraph of genome, with id colomn matching the matrix
#' @region: genome region defined above
#' @dist: distance between bins to paint, -1 means use the maximum
#' @type: ["up", "bottom", "full"], parts to paint
#' @by: axis label stride
#' @unit: axis label unit
#' @transposed: whether to transpose the original matrix (eg. up to bottom)
#' @legend_length: unit is cm
#' @line_size: size of the lines dividing chromosomes
#' @collision_mininal_distance: ticks with distance to the next tick closer
#' than this threshold will be omitted.
ggtrapezoid<-function(mat,
                       bed,
                       region=":",
                       dist=-1,
                       type="full",
                       by=5,
                       unit=1000000,
                       transposed=F,
                       legend_length=2,
                       line_size=1,
                       collision_mininal_distance=10
                       )
{
    # input transform
    colnames(bed) = c("chrx","startx","endx","idx")
    by=by*unit
    tmp = region2id(region,bed)
    Chr = tmp[3]
    start = as.integer(tmp[1])
    end = as.integer(tmp[2])
    # constants
    binSize=bed[1,endx]-bed[1,startx]
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
    datapoly[,val:=log2(val+1)]
    
    # merge to create plot data
    datapoly=datapoly[abs(y)<=dist]
    if(type=="up") datapoly=datapoly[y>=0]
    else if(type=="bottom") datapoly=datapoly[y<=0]
    
    # select 5Mb(by*unit) ticks from bed data
    xmax=2*diag_len
    tmp=bed[startx%%by==0 & idx>=start & idx<=end]
    # also add the left and right most tick
    idmin=start
    idmax=end
    tmp=rbind(bed[idx==idmin],tmp)
    tmp=rbind(tmp,bed[idx==idmax])
    # modification to avoid ugly overlapping labels
    tmp$nxt=tmp$idx[c(2:length(tmp$idx),1)]
    tmp[,dis:=nxt-idx]
    tmp=tmp[dis<0 | dis>collision_mininal_distance]
    # transformation from idx to real x coordinates.
    myBreaks=2*(tmp$idx-start)
    myLabels=round(tmp$startx/unit,digits=1)
    myChrom=data.table(breaks=myBreaks,labels=myLabels)[labels==0]
    # ggplot2
    p=ggplot(datapoly,aes(x,y))+
        geom_polygon(aes(fill=val,colour=val,group=id))+
        scale_x_continuous(position="top",
                           limits=c(0,xmax),
                           expand=c(0,0),
                           breaks=myBreaks,
                           labels=myLabels
                           )+
        
        xlab(Chr)+
        ylim(min(datapoly$y),max(datapoly$y))+
        scale_fill_gradientn(colours=c("black","yellow","orange","red","darkred"))+
        coord_fixed(ratio=1)+
        theme(panel.grid=element_blank(),
              panel.background=element_blank(),
              axis.text.y=element_blank(),
              axis.title.y=element_blank(),
              axis.ticks.y=element_blank(),
              axis.line.x = element_line(size=0.7),
              axis.ticks.length=unit(.25,'cm'),
              legend.key.height=unit(legend_length,"cm"),
              legend.title=element_blank()
              )

    # add chromosome border, if any.
    if(nrow(myChrom)!=0){
        for(i in 1:nrow(myChrom)){
            p=p+geom_abline(color="white",size=line_size,slope=1,intercept=-myChrom[i,breaks])
            p=p+geom_abline(color="white",size=line_size,slope=-1,intercept=myChrom[i,breaks])
        }
    }
    
    return(p)
}

#' paint the chromosomes
#' @bed the bedGraph file.
#' @val string of colname contains value to be painted.
ggchrom<-function(bed,val,font_size=18,border_color="white",low="gray",high="black"){
    p <- ggplot(bed,
                aes(xmin=start,xmax=end,ymin=0,ymax=1))+
         geom_rect(aes_string(fill=val))+
         facet_grid(.~chr,scales="free",space="free")+
         scale_x_continuous(expand=c(0,0))+
         theme_void()+
         theme(strip.background=element_blank(),
               strip.text.x=element_text(size=font_size),
               panel.spacing=unit(0,"cm"),
               legend.position=0,
               panel.border=element_rect(size=1,colour=border_color))+
         scale_fill_gradient(low=low,high=high)
    return(p)
}
