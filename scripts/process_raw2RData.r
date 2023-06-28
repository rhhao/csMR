#########################################################################
# File Name: process_raw2RData.r
# Author: Ruo-Han Hao and Feng Jiang
# Last Modified: 2023-01-08
# Description: convert maf and ld files to RData. Note：set sdY = 1 whem processing eQTL data.
# Usage: 
#########################################################################
library(getopt)
command<-matrix(c(
		  'help', 'h', 0,'logical', 'help.',
		  'maf', 'm', 1, 'character', 'maf file.',
		  "type",'t',1,'character','quant or cc phenotype data.',
		  "sdY",'s',2,'logical',"whether to use sdY=1."
		  ),byrow=T,ncol=5)
args<-getopt(command)
if(!is.null(args$help) ){
	help=getopt(command, usage = T)
	cat(paste(help,"convert maf and ld files to RData.\n"))
}
if(is.null(args$sdY)){
	args$sdY=FALSE
}
#library(coloc)
get_input<-function(maf_file,ld_file,bim_file,type){ 
	data=read.table(maf_file,sep="\t",header=T,row.names=1)
	bim=read.table(bim_file,sep="\t",header=F)
	data=data[as.character(bim[,2]),] 
	N=median(data[which(!is.na(data$N)),"N"]) 
	D<-list(beta=data$BETA,varbeta=data$SE^2,N=N,type=type,MAF=data$MAF,snp=rownames(data))
	#ld
	ncount=length(as.character(bim[,2]))
	ld_data<-matrix(readBin(ld_file, what="numeric", n=ncount^2, size=4),nrow = ncount,ncol = ncount)
	ld_data[is.na(ld_data)]<-0 
	rownames(ld_data)<-D$snp
	colnames(ld_data)<-D$snp
	D$LD<-ld_data
	if(args$sdY){
		D$sdY=1
	}
	return(D)
}

if(!is.null(args$maf)){
	print(paste("converting maf to RData ------------ ",args$maf," start",sep=""))
	maf=args$maf
	ld=paste(gsub("\\.maf$","",maf),".snplist.ld.bin",sep="")
	bim=paste(gsub("\\.maf$","",maf),".snplist.bim",sep="")
	D=get_input(maf,ld,bim,args$type)
	save(D,file=paste(maf,".raw.RData",sep=""))
	#check alignment
	#pdf(file=paste(maf,".check_alignment.pdf",sep=""))
	#par(mfrow=c(1,1))
	#check_alignment(D)
	#dev.off()
	print(paste("converting maf to RData ------------ ",args$maf," done",sep=""))
}
