# -*- coding: utf-8 -*-
#########################################################################
# File Name: summarize.py
# Author: Ruo-Han Hao and Feng Jiang
# Last Modified: 2023-01-08
# Description: summarize coloc results
# Usage: python summarize.py -c [coverage] -id [gwas_id] -o [outdir] -co [coloc_cutoff]
#########################################################################
import os, sys, re
import argparse

parser = argparse.ArgumentParser(description = "Summarize susie colocalization results.")
parser.add_argument('-c', '--coverage', help="Coverage used in susie finemapping.", type=str, required=True)
parser.add_argument('-id', '--gwas_id', help="GWAS ID.", type=str, required=True)
parser.add_argument('-o', '--outdir', help="The output directory.", type=str, required=True)
parser.add_argument('-co', '--coloc_cutoff', help="The coloc cutoff used to select colocalizated SNPs.", type=float, required=True)
args = parser.parse_args()
coverage_list = args.coverage.split(",")
gwas_id = args.gwas_id
outdir = args.outdir
coloc_cutoff = args.coloc_cutoff

def task(coloc_dir,coverage_list): 
    gwas_id, cell_type = re.sub("^coloc\.|\.eQTL_split\.pairs\.coloc_results$", "", coloc_dir).split(".GWAS_split.")
    fileprefix_set = set([".maf.".join(i.split(".maf.")[:2])+".maf" for i in os.listdir(coloc_dir) if re.search("\.susie\.coloc_result$",i)])
    d_res = {}
    for coverage in coverage_list:
        res=[]
        for fileprefix in fileprefix_set:
            gene_info = fileprefix.split(".maf.")[0]
            eqtl_range, ensg, cell = gene_info.split("--")
            gwas_info = re.sub("\.maf$","",fileprefix.split(".maf.")[1])
            filename=fileprefix+"."+coverage+".susie.coloc_result"
            with open(coloc_dir+"/"+filename) as f1:
                next(f1)
                for i in f1:
                    i = i.strip().split("\t")
                    i[3:8] = ["%.4g" % float(i) for i in i[3:8]]
                    res.append([gwas_id, gwas_info, eqtl_range, ensg, cell]+i)
        d_res[coverage]=res #0.9:[[line],[line],[line],]
    return d_res

def clean_up(precomp_dir, gwas_id):
    for file in os.listdir(precomp_dir):
        if re.search("eqtl_gwas_pair_file|eqtl_maf_finemap.list|eqtl_maf.list|gwas_maf_finemap.list|^coloc\.%s\.GWAS_split\..*\.eQTL_split.pairs$"%(gwas_id), file):
            os.system("rm %s"%(os.path.join(precomp_dir, file)))

##################################
precomp_dir = os.path.join(outdir, "precomputation")
outdir = os.path.join(outdir, "COLOC")
pattern=re.compile('^coloc\.'+gwas_id+'\.GWAS_split\..+\.eQTL_split\.pairs\.coloc_results$')
coloc_dir_list = [os.path.join(outdir, coloc_dir) for coloc_dir in os.listdir(outdir) if re.search(pattern, coloc_dir)]
d_res_all={}
for coloc_dir in coloc_dir_list:
    d_res = task(coloc_dir, coverage_list)
    for coverage in d_res:
        d_res_all[coverage]=d_res_all.get(coverage,[])+d_res[coverage]

for coverage in d_res_all:
    coloc_res = os.path.join(outdir,"Susie_coloc."+gwas_id+".results.coverage."+coverage)
    with open(coloc_res, 'w') as ff1:
            ff1.write("gwas_name\tgwas_info\teqtl_range\tensg\tcell\tnsnps\thit_eqtl\thit_gwas\tPP.H0.abf\tPP.H1.abf\tPP.H2.abf\tPP.H3.abf\tPP.H4.abf\tidx1\tidx2\n")
            for res in d_res_all[coverage]:
                ff1.write("\t".join([str(j) for j in res])+"\n")
    os.system("awk -F'\t' 'strtonum($13) >= %(coloc_cutoff)f' %(coloc_res)s > %(coloc_filter)s"%{"coloc_cutoff":coloc_cutoff, "coloc_res":coloc_res , "coloc_filter":coloc_res+".PPH4."+str(coloc_cutoff)})

clean_up(precomp_dir, gwas_id)
