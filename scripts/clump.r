library(dplyr)
library(data.table)
library(getopt)
command<-matrix(c(
                  'help', 'h', 0,'logical', 'help.',
                  'exp_ma_file', 'e', 1, 'character', 'input GWAS .ma file.',
                  "snp_file",'s',1,'character','colocalized SNP list file.',
                  "ref_genotype",'r',1,'character',"provide path to reference genotype.",
		  'duplication_path','d',1,'character','provide path to duplicated snp list. If not available input "None".'
                  ),byrow=T,ncol=5)

args<-getopt(command)
if(!is.null(args$help) ){
        help=getopt(command, usage = T)
        cat(paste(help,"do plink clumping from input IVs.\n"))
}

Clump <- function(exp_ma_file,ref_path,duplication_path,clump_out){
	if(duplication_path != "None"){
		for (n in 1:22){
			dup_snp_file = list.files(pattern = paste0("chr",n,"(\\D+|$)"), path = duplication_path, ignore.case=T)
			plink_file = gsub("\\.bim","",list.files(pattern = paste0("chr",n,"\\D+.*bim$"), path = ref_path, ignore.case=T))
			system(paste0('plink --silent --bfile ', ref_path, '/', plink_file,' --clump ',exp_ma_file,' --clump-field P --clump-kb 10000 --clump-p1 5e-8 --clump-p2 0.01 --clump-r2 0.001 --out ',clump_out, '.chr', n,' --exclude ', duplication_path,'/',dup_snp_file))
		}
	}else{
		for (n in 1:22){
			plink_file = gsub("\\.bim","",list.files(pattern = paste0("chr",n,"\\D+.*bim$"), path = ref_path, ignore.case=T))
			system(paste0('plink --silent --bfile ', ref_path, '/', plink_file,' --clump ',exp_ma_file,' --clump-field P --clump-kb 10000 --clump-p1 5e-8 --clump-p2 0.01 --clump-r2 0.001 --out ',clump_out, '.chr', n))
		}
	}
}

#cell <- strsplit(basename(snp_file),"PPH4\\.0\\.[1-9][0-9]*\\.")[[1]][2]
#cell <- gsub("\\.IV","",cell)

exp_ma_file = args$exp_ma_file
if (file.exists(paste0(getwd(),"/",exp_ma_file))){
	exp_ma_file = paste0(getwd(),"/",exp_ma_file)
}
snp_file = args$snp_file
ref_path = args$ref_genotype
duplication_path = args$duplication_path
outdir = paste0(dirname(snp_file),'/clump')
dir.create(outdir)
exp_ma_data <- fread(file = exp_ma_file,header = T,data.table = F)
if (file.info(snp_file)$size == 0){
	file.create(paste0(outdir,'/',gsub("\\.IV", ".clumped", basename(snp_file))))
}else{
	snp_list <- read.table(file = snp_file,header = F)[,1]
	match <- match(snp_list, exp_ma_data$SNP) %>% na.omit
	exp_ma_data <- exp_ma_data[match,]
	write.table(exp_ma_data,file = paste0(outdir,'/',gsub("\\.IV", ".ma", basename(snp_file))),row.names = F,col.names = T,sep = "\t",quote = F)
	Clump(exp_ma_file = paste0(outdir,'/',gsub("\\.IV", ".ma", basename(snp_file))), ref_path = ref_path, duplication_path = duplication_path, clump_out = paste0(dirname(snp_file),'/clump/',gsub("\\.IV", ".clump", basename(snp_file))))
	df = data.frame()
	files = list.files(outdir, pattern=paste0("^",gsub("\\.IV", "", basename(snp_file)),"\\.clump.chr.+\\.clumped$"))
	for (file in files){
		data = read.table(paste0(outdir,"/",file), header = T)
		df = rbind(df,data)
	}
	write.table(df, file = paste0(outdir,'/',gsub("\\.IV", ".clumped", basename(snp_file))), row.names = F,col.names = T,sep = "\t",quote = F)

	for (file in list.files(outdir, pattern=paste0("^",gsub("\\.IV", "", basename(snp_file)),"\\.clump\\.chr[1-9][0-9]*.*|","^",gsub("\\.IV", "", basename(snp_file)),"\\.ma$"))){
		file.remove(paste0(outdir, "/", file))
	}
}

