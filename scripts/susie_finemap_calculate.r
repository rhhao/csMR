#########################################################################
# File Name: susie_finemap_calculate.r
# Author: Ruo-Han Hao and Feng Jiang
# Last Modified: 2023-01-08
# Description: perform finemapping using susieR
# Usage: 
#########################################################################
library(coloc)
library(getopt)
command<-matrix(c(
		  'help', 'h', 0,'logical', 'help',
		  'rdata', 'r', 1, 'character', 'RData file.',
		  'type','t',1,'character','eqtl or gwas data.',
		  'gwas_type','g','1','character','GWAS data type. cc for case-control study, quant for quantitative traits',
		  'coverage',"c",1,"character",'coverage of finemapping(credible set),can split by comma'
		  ),byrow=T,ncol=5)
args<-getopt(command)
if(!is.null(args$help) ){
	help=getopt(command, usage = T)
	cat(paste(help,"perform finemapping using susieR.\n"))
}

maxit=10000

#
susie_finemap<-function(D,prior_variance,coverage){
	try_susie<-try(S <- runsusie(D,maxit=maxit,coverage=coverage, repeat_until_convergence=F))
	if('try-error' %in% class(try_susie)){ 
		try_susie2<-try(S <- runsusie(D,maxit=maxit,estimate_prior_variance=FALSE,coverage=coverage, repeat_until_convergence=F))
		if('try-error' %in% class(try_susie2)){
			try_susie3<-try(S <- runsusie(D,maxit=maxit,estimate_prior_variance=FALSE,prior_variance= prior_variance,coverage=coverage, repeat_until_convergence=F))
			if('try-error' %in% class(try_susie3)){
				S <- "NA"
				print("All susie failed")
				q()
			}
			else{
				print("Try susie3 succeed")
			}
		}
		else{
			print("Try susie2 succeed")
		}
	}
	else{
		print("Try susie1 succeed")
	}
	return(S)
}

task<-function(coverage){
	if(args$type=="eqtl"){
		S <- susie_finemap(D,(0.15/1)^2,coverage)
	}
	else if(args$type=="gwas"){
		if(args$gwas_type=="cc"){
			S <- susie_finemap(D,0.2^2,coverage)
		}else if(args$gwas_type=="quant"){
			S <- susie_finemap(D,(0.15/1)^2,coverage)
		}else{
			print("please select gwas data type \"cc\" or \"quant\"")
			q()
		}

	}
	else{
		print("please select type \"eqtl\" or \"gwas\"")
		q()
	}
	if(!identical(S,"NA")){
		maf=gsub("\\.raw\\.RData$","",args$rdata)
		susie_name=paste(maf,".",coverage,".susie.RData",sep="")
		save(S,file=susie_name)
		print(paste("finemapping is succeed and done! ",susie_name,sep=""))
	}
	else{
		print(paste("finemapping is failed! ",susie_name,sep=""))
	}

}

if(!is.null(args$rdata)){
	print(paste("susie finemapping ------------ ",args$rdata,sep=""))
	coverage_list=strsplit(args$coverage,",")[[1]]
	load(args$rdata)
	for(i in coverage_list){
		task(i)
	}
	}
