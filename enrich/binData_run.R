#! /usr/bin/env Rscript


Args<-commandArgs(trailingOnly = T)

if (length(Args)==0) {
    stop("At least 1 argument must be supplied: Rscriipt script.R dssInput.gz", call.=FALSE)
}

## Import necessary packages ----
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(Rcpp))

sourceCpp("bin100_RcppVer.cc")
source("calculate_binom.R")

fread(paste0("zcat < ",Args[1]))->data
colnames(data)=c("chr","pos","N","X")
calculateBinomTest(data)->data
calculateBinMethylation(data,100)->bin100
calculateBinMethylation(data,2000)->bin2k
fwrite(bin100,str_c(Args[1],".bin100.csv"),sep=",")
fwrite(bin2k,str_c(Args[1],".bin2k.csv"),sep=",")


