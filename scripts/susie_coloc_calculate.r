#########################################################################
# File Name: susie_coloc_calculate.r
# Author: Ruo-Han Hao and Feng Jiang
# Last Modified: 2023-01-08
# Description: perform colocalization by susie
# Usage: 
#########################################################################
library(coloc)
library(getopt)
command<-matrix(c(
		  'help', 'h', 0,'logical', 'help',
		  'eqtl_maf', 'e', 1, 'character', 'eqtl maf file',
		  'gwas_maf', 'g', 1, 'character', 'gwas maf file',
		  'coverage', 'c', 2, 'character', 'susie finemapping coverage',
		  'out_dir', 'o', 1, 'character', 'output directory'),byrow=T,ncol=5)
args<-getopt(command)
if(!is.null(args$help) ){
	help=getopt(command, usage = T)
	cat(paste(help,"perform colocalization by susie\n"))
}
#if coverage is not defined use default coverage 0.9
if(is.null(args$coverage)){
	args$coverage="0.9"
}

if(!is.null(args$eqtl_maf)){
	print(paste("process susie_coloc ------------ ",args$eqtl_maf," ------- ",args$gwas_maf,sep=""))
	out_dir=args$out_dir
	eqtl_maf=args$eqtl_maf
	gwas_maf=args$gwas_maf
	coverage=args$coverage
	#eqtl_raw_file=paste(eqtl_maf,".raw.RData",sep="")
	#gwas_raw_file=paste(gwas_maf,".raw.RData",sep="")
	eqtl_susie_file=paste(eqtl_maf,".",coverage,".susie.RData",sep="")
	gwas_susie_file=paste(gwas_maf,".",coverage,".susie.RData",sep="")
	#load data
	#load(eqtl_raw_file)
	#D_eqtl<-D
	#load(gwas_raw_file)
	#D_gwas<-D
	load(eqtl_susie_file)
	S_eqtl<-S
	load(gwas_susie_file)
	S_gwas<-S
	#n_intersect=length(intersect(D_eqtl$snp,D_gwas$snp))
	n_intersect=length(intersect(names(S_eqtl$pip),names(S_gwas$pip)))
	if(n_intersect<2){
		print(paste(eqtl_maf,gwas_maf,"JFF coloc only intersect 0 or 1 SNPs",sep="\t"))
		q()
	}

	# susie coloc
	susie_res <- coloc.susie(dataset1=S_eqtl, dataset2=S_gwas,trim_by_posterior=FALSE)
	write.table(susie_res$summary,file=paste(out_dir,"/",basename(gwas_maf),".",basename(eqtl_maf),".",coverage,".susie.coloc_result",sep=""),quote = FALSE, sep = "\t",row.names=FALSE)
	#save(susie_res,file=paste(out_dir,"/",basename(gwas_maf),".",basename(eqtl_maf),".",coverage,".susie.coloc_result.RData",sep=""))
	#sensitivity test
	#D_eqtl$position=seq(length(D_eqtl$snp))
	#D_gwas$position=seq(length(D_gwas$snp))
	
	#n=dim(susie_res$summary)[1]
	#pdf(file=paste(out_dir,"/",basename(gwas_maf),".",basename(eqtl_maf),".",coverage,".susie.sensitivity.pdf",sep=""))
	#par(mfrow=c(2,2))
	#for(i in seq(n)){
	#		sensitivity(susie_res,"H4 > 0.8",row=i,dataset1=D_eqtl,dataset2=D_gwas,preserve.par=TRUE)
	#}
	#dev.off()
}

