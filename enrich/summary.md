## Boxplot

1. batch_enrich.sh

   ```bash
   #! /bin/bash
   for f in *.bedGraph
   do
       Rscript enrich.R gr_TE_chr.rds interestedTEcls.rds $f ${f%%.bedGraph}
   done
   ```

2. enrich.R

   ```R
   #! /usr/bin/env Rscript
   Args<-commandArgs(trailingOnly = T)
   
   # ---------------
   # enrich.R: plot boxplot for each histone modification and output plot data.
   # totalRegion.rds: data.table of all the regions. (here are total TEs)
   # interestedRegino.rds: GenomicRanges object of interested regions from different classes.
   # 	(here are 5 cls of Differentially expressed-TEs)
   # binData.bedGraph: Histone bedGraph file.
   # outputName: prefix for output files.
   # <!> the input's chromosome format (chrM/Mt/ChrM) must be the same.
   # ---------------
   
   if (length(Args)!=4) {
       stop("Usage:Rscript enrich.R totalRegion.rds interestedRegion.rds binData.bedGraph outputName", call.=F)
   }
   
   # Import necessary packages ----
   suppressPackageStartupMessages(library(GenomicRanges))
   suppressPackageStartupMessages(library(data.table))
   suppressPackageStartupMessages(library(stringr))
   suppressPackageStartupMessages(library(Rcpp))
   suppressPackageStartupMessages(library(ggplot2))
   suppressPackageStartupMessages(library(ggsignif))
   
   sourceCpp("sumBinOverRegion.cc")
   
   readRDS(Args[1])->regions
   readRDS(Args[2])->interestedRegions
   fread(Args[3])->binData
   outputName=Args[4]
   
   colnames(binData)=c("seqnames","start","end","val")
   as.data.table(regions)->regions
   as.data.table(interestedRegions)->interestedRegions
   
   regions$cls="genome"
   # accumulate bin value in specific regions.
   regions$val=sumBinOverRegion(regions,binData)
   # normalize with TE(region) width.
   regions[,val:=val/width]
   
   interestedRegions$val=sumBinOverRegion(interestedRegions,binData)
   interestedRegions[,val:=val/width]
   
   # rbind and extract pdata
   rbind(regions,interestedRegions)[,c("cls","val")]->pdata
   
   ggplot(pdata,aes(as.factor(cls),val,fill=as.factor(cls)))+
       geom_boxplot()+
       geom_jitter(data=pdata[cls!="genome"],position=position_jitter(0.2),alpha=0.5,size=3)+
       scale_y_log10()+
       xlab("class")+
       ylab("density")+
       ggtitle(str_c(str_replace(outputName,"(.*)-.*","\\1")," Enrichment"))+
       theme_bw()+
       theme(legend.position="none",
             text=element_text(size=16),
             plot.title=element_text(hjust=0.5))+
       geom_signif(comparisons=list(c("genome","1"),
                                    c("genome","2"),
                                    c("genome","3"),
                                    c("genome","4"),
                                    c("genome","5")),
                   map_signif_level=T,step_increase=0.1)
   
   # output plot data for heatmap
   fwrite(pdata,str_c(outputName,"_pdata.csv",sep=","))
   fwrite(interestedRegions,str_c(outputName,".csv"),sep=",")
   ggsave(str_c(outputName,".png"))
   
   ```

3. sumBinOverRegions.cc

   ```c++
   #include <Rcpp.h>
   #include <string>
   #include <cstring>
   #include <vector>
   using namespace Rcpp;
   //number of regions to process at most
   const int region_max = 40000;
   /*
   sumBinOverRegion: 
   regions: targeted regions to calculate accumulation.
   bins: binData across the genome.
   */
   // [[Rcpp::export]]
   NumericVector sumBinOverRegion(DataFrame regions, DataFrame bins, bool verbose=0){
       //convert data.frame to vector
       CharacterVector chr_region=regions["seqnames"];
       IntegerVector start_region=regions["start"];
       IntegerVector end_region=regions["end"];
   
       CharacterVector chr_bin=bins["seqnames"];
       IntegerVector start_bin=bins["start"];
       IntegerVector end_bin=bins["end"];
       NumericVector val=bins["val"];
   	
       //array of accumulated values (result)
       double sval[region_max];
       memset(sval,0,sizeof(sval));
   
       int len = chr_bin.size();
       int rlen = chr_region.size();
       if(verbose) Rcout<<"len of binData:"<<len<<"\n";
   
       int r=0;  //id of current region
       //travel through the genome(binData)
       for(int i=0;i<len;i++){
           if(verbose && i%100000==0) Rcout<<i<<" "<<r<<"\n";
           //for variable binSize which may cover a whole region
           if(chr_region[r]==chr_bin[i] &&
              start_bin[i]<=start_region[r] &&
              end_bin[i]>=end_region[r])
           {
               if(verbose) Rcout<<"region: "<<r<<" is too short"<<"\n";
               sval[r]+=val[i];
               r++;
               if(r==rlen) break;
           }
           //normal accumulation
           else if(chr_region[r]==chr_bin[i] &&
                   end_bin[i] >= start_region[r] &&
                   start_bin[i] <= end_region[r])
           {
               sval[r]+=val[i];                    
               if(end_bin[i] > end_region[r]){
                   if(verbose) 
                       Rcout<<"finished regions to:"<<r<<
                       "("<<chr_region[r]<<": "<<start_region[r]<<","<<end_region[r]<<")"<<
                       ", sum is "<<sval[r]<<"\n";
                   r++;
                   if(r==rlen) break;
                   //for neighboured regions and regions inside regions: slide back.
                   while(end_bin[i]>=start_region[r]){
                       i--;
                   }
               }
           }
           //for some possible gaps in binData
           while(chr_region[r]==chr_bin[i] && start_bin[i]>end_region[r]){
               if(verbose) Rcout<<"region: "<<r<<" is jumped over\n";
               r++;
               if(r==rlen) break;
           }
           if(r==rlen) break;
       }
   
       if(verbose)Rcout<<"processed regions: "<< r <<"\n";
       return wrap(std::vector<double>(sval,sval+r));
   }
   
   ```

   

## Heatmap

```R
#input data
dir()[str_detect(dir(),"pdata")]->files
fread("H2A.W-SRX352374,_pdata.csv")->pdata

# cbind and change colname for each histone modification
colnames(pdata)[ncol(pdata)]="H2A.W"
for(f in files){
    cbind(pdata,fread(f)[,val])->pdata;
    colnames(pdata)[ncol(pdata)]=str_replace(f,"(.*)-.*","\\1")
}

#example
if(F){
pdata
          cls      H2A.W      H2A.X      H2A.Z         H3.1         H3.3
    1: genome 0.04262500 0.14912500 0.24875000 0.1203750000 0.0672500000
    2: genome 0.06535433 0.10566929 0.04897638 0.0000000000 0.0000000000
    3: genome 0.12087849 0.06872173 0.01804840 0.0000000000 0.0008048396
    4: genome 0.25592949 0.06272436 0.01798077 0.0000000000 0.0000000000
    5: genome 0.12356271 0.07481687 0.01049945 0.0087236404 0.0227635960
   ---
31433:      1 0.17971529 0.07903198 0.01959702 0.0072798949 0.0049277267
31434:      1 0.03056196 0.11494236 0.10017291 0.0006195965 0.0014697406
31435:      1 0.02544776 0.18089552 0.06589552 0.0097014925 0.0000000000
31436:      4 0.41818644 0.15698305 0.03616949 0.0062796610 0.0067118644
31437:      1 0.09875000 0.12704545 0.10915909 0.0218636364 0.3964318182
}

#calculate mean of each class
apply(pdata[cls=="genome",-1],2,mean)->pdata2
for(t in 1:5){
    rbind(pdata2,apply(pdata[cls==t,-1],2,mean))->pdata2
}
rownames(pdata2)=c("genome",1:5)

#apply some transformations
t(pdata2)->pdata2t
as.data.table(pdata2t,keep.rownames=T)->pdata2t

#apply log2(histone/genome)
log2(pdata2t[,2:7]/pdata2t$genome)->pdt_divlog2
pdt_divlog2[,rn:=pdata2t$rn]
pdt_divlog2[,genome:=NULL]

#example
if(F){
 pdt_divlog2
        class1      class2     class3       class4      class5       rn
 1:  0.5101522  0.58538496  0.2967470  0.302612896 -0.06148577    H2A.W
 2:  0.1365734  0.23043011  0.1959248  0.401607484  0.14911672    H2A.X
 3:  0.7579350 -0.19767231  1.0240661 -0.006585763  0.74382338    H2A.Z
 4:  0.5787154  0.55145716  0.2477744 -0.371900540 -2.04371814     H3.1
 5:  0.6071332  0.08530672 -0.2110158 -0.310798840 -1.22241817     H3.3
 6:  0.1053927 -0.25920153 -0.1850042 -0.159205207  0.05385622  H3K14ac
 7:  0.4169843  0.04757198 -0.2561090  0.007392206  0.33654714  H3K23ac
 8:  0.6543722  0.22877841  0.2842896  0.020284815  0.71566648  H3K27ac
 9:  1.4996341  1.71939857  0.7432015  1.779393061 -0.10202195 H3K27me1
10: -0.1263378 -1.29681958 -0.5863218 -0.120452033 -0.19439334 H3K27me3
11:  0.6262454 -0.18543351 -0.3516174 -0.486306268  0.26257599  H3K36ac
12:  0.2392424  0.23563357  1.3865974 -0.435059504 -0.06320146 H3K36me3
13:  1.5804816  1.72442521  2.8293885  1.346526153  0.04926893  H3K4me1
14:  1.4844249  0.12711697  1.3281237  0.055389048  1.04872611  H3K4me2
15:  0.4992290 -0.01659348  0.2353465 -0.092735039  0.60648168  H3K4me3
16:  0.2666571 -0.70663652 -0.4636182  0.211508066  0.82230263  H3K56ac
17: -0.1325705 -0.54706379 -0.7613678 -0.968375205 -0.03932836   H3K9ac
18:  0.8664797  0.88931363  0.4830381  0.968403072  0.05438681  H3K9me1
19:  0.8729202  1.08786180  0.7124731  1.285640470  0.40152427  H3K9me2
20:  0.5369150 -0.23046708  0.3365756 -0.473974292  0.91560190  H4K16ac
}

# calculate p_value of wilcox.test between genome and other 5 types
pdt_pvalue=pdt_divlog2
for(t in 1:5){
    for(his in pdt_pvalue$rn){
    	wilcox.test(pdata[cls=="genome"][[his]],pdata[cls==t][[his]])$p.value->pdt_pvalue[rn==his,t]
    }
}

# assign 0 to not significant pairs. 
pdt_divlog2_pcut=pdt_divlog2
for(i in 1:nrow(pdt_divlog2_pcut)){
    for(j in 1:(ncol(pdt_divlog2_pcut)-1)){
        if(pdt_pvalue[i][[j]]>0.05) pdt_divlog2_pcut[i][[j]]=0
    }
}

# some transformations for plotting
as.data.frame(pdt_divlog2_pcut[,1:5])->pdf_dvlg2_pct5
rownames(pdf_dvlg2_pct5)=pdt_divlog2$rn

# set color of 0 point to white
mycolor=colorRampPalette(c("blue","white","red"))(50)
mybreaks=c(seq(min(pdf_dvlg2_pct5),0,length.out=26),
          seq(max(pdf_dvlg2_pct5)/50,max(pdf_dvlg2_pct5),length.out=25))

#example
if(F){
pdf_dvlg2_pct5
             class1      class2     class3      class4      class5
H2A.W     0.5101522  0.58538496  0.0000000  0.30261290  0.00000000
H2A.X     0.1365734  0.23043011  0.0000000  0.40160748  0.14911672
H2A.Z     0.0000000  0.00000000  0.0000000  0.00000000  0.00000000
H3.1      0.5787154  0.55145716  0.2477744 -0.37190054 -2.04371814
H3.3      0.6071332  0.08530672 -0.2110158 -0.31079884 -1.22241817
H3K14ac   0.0000000 -0.25920153  0.0000000  0.00000000  0.00000000
H3K23ac   0.4169843  0.00000000  0.0000000  0.00000000  0.00000000
H3K27ac   0.6543722  0.22877841  0.0000000  0.00000000  0.71566648
H3K27me1  1.4996341  1.71939857  0.0000000  1.77939306 -0.10202195
H3K27me3  0.0000000  0.00000000  0.0000000  0.00000000  0.00000000
H3K36ac   0.6262454 -0.18543351  0.0000000  0.00000000  0.00000000
H3K36me3  0.2392424  0.00000000  0.0000000  0.00000000  0.00000000
H3K4me1   1.5804816  1.72442521  0.0000000  1.34652615  0.04926893
H3K4me2   1.4844249  0.00000000  0.0000000  0.00000000  1.04872611
H3K4me3   0.4992290 -0.01659348  0.2353465 -0.09273504  0.60648168
H3K56ac   0.2666571  0.00000000  0.0000000  0.00000000  0.00000000
H3K9ac   -0.1325705 -0.54706379  0.0000000 -0.96837521  0.00000000
H3K9me1   0.8664797  0.88931363  0.0000000  0.96840307  0.00000000
H3K9me2   0.8729202  1.08786180  0.0000000  1.28564047  0.00000000
H4K16ac   0.0000000  0.00000000  0.0000000 -0.47397429  0.00000000
}

pheatmap(pdf_dvlg2_pct5,cluster_rows=F,cluster_cols=F,border_color="white",color=mycolor,breaks=mybreaks)

dev.print(pdf,"enrichment.pdf")
```





