library(ggplot2)
library(data.table)

ggtrapezoid<-function(mat, start, end, dist, flip=F, type="full"){
    #flip the matrix, considering the mat is in the format of the previously
    #viewed heatmap. If my consideration is somewhat wrong, just comment it
    #and try again.
    if(flip) mat=mat[c(nrow(mat):1),]
    # constants
    diag_len = end - start + 1
    total_num = diag_len*diag_len
    # error input 
    if(nrow(mat) != ncol(mat) |
       diag_len > nrow(mat) |
       dist > diag_len )
    {
        return("Input data error!")
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
    mat2=mat[start:end,start:end]
    # for debug, print what we plot (t(mat2) clockwise 45 degree --> plot)
    # print(t(mat2))
    # flatten the matrix for plotting
    values=data.table(id=1:total_num,val=as.vector(as.matrix(t(mat2))))
    
    # merge to create plot data
    values[diamonds_vertices,on=("id")]->datapoly
    datapoly=datapoly[abs(y)<=dist]
    if(type=="up") datapoly=datapoly[y>=0]
    else if(type=="bottum") datapoly=datapoly[y<=0]
    
    # ggplot2
    p=ggplot(datapoly,aes(x,y))+
        geom_polygon(aes(fill=val,group=id))+
        xlim(0,2*diag_len-1)+
        ylim(min(datapoly$y),max(datapoly$y))+
        scale_fill_gradientn(colours=c("white","red"))+
        coord_fixed(ratio=1)

    return(p)
}
