library(dplyr)
library(data.table)
library(getopt)
library(phenoscanner)

command<-matrix(c(
                  'help', 'h', 0,'logical', 'help.',
                  'coloc_res', 'r', 1, 'character', 'susie coloc result.',
                  'exp_ma_file', 'e', 1, 'character', 'input GWAS .ma file.',
                  'cells','c',1,'character','all tested cells.',
                  'outdir','o',1,'character',"output directory."
                  ),byrow=T,ncol=5)
args<-getopt(command)
if(!is.null(args$help) ){
        help=getopt(command, usage = T)
        cat(paste(help,"select MR IVs.\n"))
}
coloc_res = args$coloc_res
exp_ma_file = args$exp_ma_file
if (file.exists(paste0(getwd(),"/",exp_ma_file))){
	exp_ma_file = paste0(getwd(),"/",exp_ma_file)
}
cells = args$cells
outdir = args$outdir
################################################
outdir = paste0(outdir, "/MR")
if (!dir.exists(paste0(outdir))){
dir.create(paste0(outdir))
}
tryCatch({
	data <- read.table(coloc_res, header = FALSE)
}, error = function(e) {
	print("An error occurred while reading the colocalization result file. Could be no colocalized SNP under PPH4 > 0.8")
	data <- data.frame(matrix(nrow = 0, ncol = 15))
	colnames(data) <- c("gwas_name", "gwas_info", "eqtl_range", "ensg", "cell", "nsnps", "hit_eqtl", "hit_gwas", "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf", "idx1", "idx2")
}, finally = {
	colnames(data) <- c("gwas_name", "gwas_info", "eqtl_range", "ensg", "cell", "nsnps", "hit_eqtl", "hit_gwas", "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf", "idx1", "idx2")
})
confounders <- c("smoking","alcohol","education","college","university","smoke","drink","drinking")
exp_ma_data <- fread(file = exp_ma_file, header = T, data.table = F, select = c("SNP", "P"))
exp_ma_data <- exp_ma_data[na.omit(match(data[,8],exp_ma_data$SNP)),] %>% unique()
exp_ma_data <- exp_ma_data[exp_ma_data$P < 5e-8,]
if(nrow(exp_ma_data) == 0){
	confounder_snps=c()
	print("Notice: No availuable IVs!")
}else if(nrow(exp_ma_data) <= 100 & nrow(exp_ma_data) > 0){
	scan_res <- phenoscanner(snpquery = exp_ma_data$SNP, pvalue = 5e-8, r2 = 0.8)$results
	confounder_snps <- scan_res[grep(paste(confounders,collapse="|"),tolower(scan_res$trait)),"rsid"] %>% unique()
}else{
	confounder_snps=c()
	for(i in seq(floor(nrow(exp_ma_data)/100))){
		scan_res <- phenoscanner(snpquery = exp_ma_data$SNP[(100*(i-1)+1):(100*i)], pvalue = 5e-8)$results
		snps = scan_res[grep(paste(confounders,collapse="|"),tolower(scan_res$trait)),"rsid"] %>% unique()
		confounder_snps = c(confounder_snps, snps)
	}
	if(nrow(exp_ma_data)%%100 != 0){
		scan_res <- phenoscanner(snpquery = exp_ma_data$SNP[(nrow(exp_ma_data)-(nrow(exp_ma_data)%%100)+1):nrow(exp_ma_data)], pvalue = 5e-8)$results
		snps = scan_res[grep(paste(confounders,collapse="|"),tolower(scan_res$trait)),"rsid"] %>% unique()
		confounder_snps = c(confounder_snps, snps)
	}
}
if(length(confounder_snps) != 0){
	exp_ma_data <- exp_ma_data[-which(exp_ma_data$SNP %in% confounder_snps),]
}
data2 <- merge(exp_ma_data, data, by.x="SNP", by.y="hit_gwas")
for (cell in strsplit(cells, split=",")[[1]]){
	if (! cell %in% data2$cell){
		file.create(paste0(outdir, "/", basename(coloc_res), ".", cell, ".IV"))
	}else{
		snplist <- data2[data2$cell == cell, "SNP"]
		write.table(unique(snplist), file = paste0(outdir, "/", basename(coloc_res), ".", cell, ".IV"), sep="\n", quote = F, col.names = F, row.names = F)
	}
}
l = paste(outdir, list.files(outdir, pattern = paste0(basename(coloc_res),".+.IV")), sep = "/")
write.table(l, file = paste0(outdir, "/", basename(coloc_res), ".", "cells"), sep="\n", quote = F, col.names = F, row.names = F)
