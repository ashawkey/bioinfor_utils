#include <Rcpp.h>
#include <string>
#include <cstring>
#include <vector>
using namespace Rcpp;
const int region_max = 40000;
/*
 * findOverlaps and sum val in special regions.
 * 
 * */

// [[Rcpp::export]]
NumericVector sumBinOverRegion(DataFrame regions, DataFrame bins, bool verbose=0){
    CharacterVector chr_region=regions["seqnames"];
    IntegerVector start_region=regions["start"];
    IntegerVector end_region=regions["end"];

    CharacterVector chr_bin=bins["seqnames"];
    IntegerVector start_bin=bins["start"];
    IntegerVector end_bin=bins["end"];
    NumericVector val=bins["val"];

    double sval[region_max];
    memset(sval,0,sizeof(sval));

    int len = chr_bin.size();
    int rlen = chr_region.size();
    if(verbose) Rcout<<"len of binData:"<<len<<"\n";
      

    int r=0;  //id_region
    for(int i=0;i<len;i++){
        if(verbose && i%100000==0) Rcout<<i<<" "<<r<<"\n";
        //for variable binSize which may cover a whole region
        if(chr_region[r]==chr_bin[i] && start_bin[i]<=start_region[r] && end_bin[i]>=end_region[r]){
            if(verbose) Rcout<<"region: "<<r<<" is too short"<<"\n";
            sval[r]+=val[i];
            r++;
            if(r==rlen) break;
        }
        else if(chr_region[r]==chr_bin[i] && end_bin[i] >= start_region[r] && start_bin[i] <= end_region[r]){
            sval[r]+=val[i];                    
            if(end_bin[i] > end_region[r]){
                if(verbose) Rcout<<"finished regions to:"<<r<<"("<<chr_region[r]<<": "<<start_region[r]<<","<<end_region[r]<<")"<<", sum is "<<sval[r]<<"\n";
                r++;
                if(r==rlen) break;
                //for neighboured regions and regions inside regions, slide
                //back bins
                while(end_bin[i]>=start_region[r]){
                    i--;
                }
            }
        }
        //for some possible uncontinuous bin jumping...
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
