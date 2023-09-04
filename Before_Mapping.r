#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = T)

library(data.table)
library(tidyr)

tresh=as.numeric(args[2])
geno=fread(args[1],data.table=F,header=F)
row.names(geno)=geno[,1]
geno=geno[,-1]


genoComp=geno

genoComp[genoComp=="A"]=1
genoComp[genoComp=="B"]=2
genoComp[genoComp=="C"]=3
genoComp[genoComp=="D"]=4
genoComp[genoComp=="H"]=5
genoComp[genoComp=="a"]=1
genoComp[genoComp=="b"]=2
genoComp[genoComp=="c"]=3
genoComp[genoComp=="d"]=4
genoComp[genoComp=="h"]=5

genoComp=genoComp %>% unite('xxx',sep='')
fwrite(genoComp,'Converted_Haplotype',row.names=T,col.names=F,sep='\t')

genoComp=geno

genoComp[genoComp=="A"]=3
genoComp[genoComp=="B"]=-3
genoComp[genoComp=="C"]=-2
genoComp[genoComp=="D"]=2
genoComp[genoComp=="H"]=1
genoComp[genoComp=="a"]=3
genoComp[genoComp=="b"]=-3
genoComp[genoComp=="c"]=-2
genoComp[genoComp=="d"]=2
genoComp[genoComp=="h"]=1
genoComp[genoComp=="-"]=NA
genoComp=apply(genoComp,1,as.numeric)


j=0
repulsion=NULL
FinalMap=NULL

while(ncol(genoComp)>2) {
	tmp=genoComp[,1]*genoComp[,2:ncol(genoComp)]
	tmp[tmp==3]=-3
	x=apply(tmp,2,function(y) length(which(y<=-3)))
	if(length(which(x<=tresh))>0) repulsion=rbind(repulsion,cbind(colnames(genoComp)[1],as.data.frame(x[which(x<=tresh)])))

	FinalMap=c(FinalMap,colnames(genoComp)[1])
	dlt=c(1,(which(x<=tresh)+1))
	genoComp=genoComp[,-dlt]
}
MissingSNP=dim(geno)[1]-dim(repulsion)[1]-length(FinalMap)
if(MissingSNP>0) FinalMap=c(FinalMap,row.names(geno)[(nrow(geno)-MissingSNP+1):nrow(geno)])
FinalMap=geno[FinalMap,]
fwrite(FinalMap,'For_Mapping.txt',row.names=T,col.names=F,sep='\t')
fwrite(repulsion,'Repulsions.txt',row.names=T,col.names=F,sep='\t')
