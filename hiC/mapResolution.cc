/* * * * * * * * * * * * * * * * * * * * * * * * * *
 * Used for calculation of HiC resolution.
 * Please change the 'total' variable to targeted genome size.
 * (default is Arabidopsis thaliana)
 * * * * * * * * * * * * * * * * * * * * * * * * * */
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <cmath>
#define ll long long

using namespace std;

char* inputFileName;
char* outputFileName;
fstream input, output;
const int contactThreshold=1000;
const int stride=3000;
const ll total=119667750;
int binsValid,binsTotal,threshold;
int binSize=50;

struct chrpos{
    string chr;
    int pos;
    chrpos(string c, int p):chr(c),pos(p) {}
    bool operator<(const chrpos& cp2) const{
        return (chr==cp2.chr)?(pos<cp2.pos):chr<cp2.chr;
    }
};

map<chrpos,int> cov50,cov;

void calc(bool isEmpty){
    cov50.clear();
    if(isEmpty){
        string read,chr1,str1,chr2,str2,ex1,ex2,ex3,ex4;
        int pos1,pos2,frag;
        ll cnt=0;
        //int qua1,qua2,frag1,frag2;
        //while(input>>str1>>chr1>>pos1>>frag1>>str2>>chr2>>pos2>>frag2>>qua1>>ex1>>ex1>>qua2>>ex1>>ex1>>ex1>>ex1){
        while(input>>read>>chr1>>pos1>>str1>>chr2>>pos2>>str2>>frag>>ex1>>ex2>>ex3>>ex4){
            //if(qua1==0||qua2==0) continue;
            //if(frag1==frag2) continue;
            cnt++;
            if(cnt%1000000==0) cout<<"read in lines: "<< cnt <<" * 1000000"<<endl;
            chrpos cp1(chr1,pos1/binSize*binSize),cp2(chr2,pos2/binSize*binSize);
            cov50[cp1]+=1;
            cov50[cp2]+=1;
        }
    }
    else{
        string chr;
        int pos, cnt;
        while(input>>chr>>pos>>cnt){
            chrpos cp(chr,pos);
            cov50[cp]=cnt;
        }
    }
    cout<<"Read in finished"<<endl;
    binsValid=0;
    for(pair<chrpos,int> p:cov50){
        chrpos& cp=p.first;
        if(isEmpty) {
            output<<cp.chr<<" "<<cp.pos<<" "<<p.second<<endl;
        }
        if(p.second>contactThreshold){
            binsValid++;
        }
    }
    if(isEmpty) cout<<"50bp coverage written out"<<endl;
    cout<<"strided search started"<<endl;
    binsTotal=total/50;
    threshold=binsTotal*4/5;
    int low=0,high=0;
    //+1000bp search
    while(binsValid < threshold){
        low=binSize;
        binSize=binSize+stride;
        cov.clear();
        binsValid=0;
        for(pair<chrpos,int> p:cov50){
            chrpos cp(p.first.chr,p.first.pos/binSize*binSize);
            cov[cp]+=p.second;
        }
        for(pair<chrpos,int> p:cov){
            if(p.second>contactThreshold){
                binsValid++;
            }
        }
        binsTotal=total/binSize;
        threshold=binsTotal*4/5;
    }
    cout<<"binary search started"<<endl;
    high=binSize;
    cout<<"low:"<<low<<" high:"<<high<<endl;
    binSize=low+(high-low)/2;
    binSize=(binSize+49)/50*50;
    //binary search from 50+n*stride to 50+(n+1)*stride
    while(binSize<high){
        cov.clear();
        binsValid=0;
        for(pair<chrpos,int> p:cov50){
            chrpos cp(p.first.chr,p.first.pos/binSize*binSize);
            cov[cp]+=p.second;
        }
        for(pair<chrpos,int> p:cov){
            if(p.second>contactThreshold){
                binsValid++;
            }
        }
        binsTotal=total/binSize;
        threshold=binsTotal*4/5;
        if(binsValid<threshold){
            low=binSize;
            binSize=low+(high-low)/2;
            binSize=(binSize+49)/50*50;
        }
        else{
            high=binSize;
            binSize=low+(high-low)/2;
            binSize=(binSize+49)/50*50;
        }
    }
    cout<<"The map resolution is "<<binSize<<endl;
}


int main(int argc,char** argv){
    ios::sync_with_stdio(false);
    if(argc==3){
        inputFileName=argv[1];
        outputFileName=argv[2];
    }
    else if(argc==2 && argv[1]!="-h" && argv[1]!="--help"){
        inputFileName=argv[1];
    }
    else{
        cout<<"Usage: ./mapResolution Hic.input Cov50.output"<<endl;
        cout<<"Or: ./mapResolution Cov50.input"<<endl;
        cout<<"Hic.input file format:"<<endl;
        cout<<"    read name / chr_reads1 / pos_reads1 / strand_reads1 / chr_reads2 / pos_reads2 / strand_reads2 / fragment_size 4*[/ allele_specific_tag]"<<endl;
        exit(1);
    }
    input.open(inputFileName,ios::in);
    if(!input.is_open()){
        cout<<"Error: can't open input file"<<endl;
        exit(1);
    }
    if(argc==3){
        output.open(outputFileName,ios::out);
        if(!output.is_open()){
            cout<<"Error: can't open output file:"<<outputFileName<<endl;
            exit(1);
        }
        calc(true);
        output.close();
    }
    else if(argc==2){
        cout<<"Use the previously generated cov50 file"<<endl;
        calc(false);
    }
    input.close();
    return 0;
}
