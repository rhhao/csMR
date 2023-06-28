# -*- coding: utf-8 -*-
#########################################################################
# File Name: susie_finemap.py
# Author: Ruo-Han Hao and Feng Jiang
# Last Modified: 2023-01-08
# Description: perform finemapping using susieR
# Usage: python susie_finemap.py -o [outdir] -d [gwas_dict] -c [coverage] -j [threads]
#########################################################################
import os,sys
import json
import re
import argparse
from multiprocessing import Pool
parser = argparse.ArgumentParser(description="Susie finemap for all RData.")
parser.add_argument('-o', '--outdir', help="The ouput directory", type=str, required=True)
parser.add_argument('-d', '--dict', help="The dictionary of GWAS types", type=str, required=True)
parser.add_argument('-c', '--coverage', help="coverage, default 0.9, can split by comma", type=str, default=0.9, required=False)
parser.add_argument('-j', '--threads', help="Threads number, default 8", type=int, default=8, required=False)
args = parser.parse_args()
outdir = args.outdir
gwas_dict = args.dict
gwas_dict = json.loads(gwas_dict)
coverage = args.coverage
threads = args.threads

def get_maf_list(outdir, coverage):
    outdir = os.path.join(outdir, "precomputation")
    pair_list = [os.path.join(outdir, file) for file in os.listdir(outdir) if re.search('^coloc\..*\.GWAS_split\..*\.eQTL_split\.pairs$', file)]
    eqtl_maf_list=[]
    gwas_maf_list=[]
    for i in pair_list:
        with open(i) as f:
            for line in f:
                gwas_maf = line.strip().split("\t")[0]
                eqtl_maf = line.strip().split("\t")[1]
                if not os.path.exists(eqtl_maf+"."+coverage+".susie.RData"):
                    eqtl_maf_list.append(eqtl_maf)
                if not os.path.exists(gwas_maf+"."+coverage+".susie.RData"):
                    gwas_maf_list.append(gwas_maf)
    eqtl_maf_list=set(eqtl_maf_list)
    gwas_maf_list=set(gwas_maf_list)
    with open(os.path.join(outdir,"gwas_maf_finemap.list"), 'w') as fw:
        fw.write("\n".join(gwas_maf_list)+"\n")
    with open(os.path.join(outdir,"eqtl_maf_finemap.list"), 'w') as fw:
        fw.write("\n".join(eqtl_maf_list)+"\n")
    return([list(gwas_maf_list), list(eqtl_maf_list)])


def task_susie(rdata_path, data_type, gwas_type): 
    os.system("Rscript scripts/susie_finemap_calculate.r -r %(rdata_path)s -t %(data_type)s -c %(coverage)s -g %(gwas_type)s" % {"rdata_path": rdata_path, "data_type": data_type,"coverage":coverage, "gwas_type":gwas_type})
    if os.path.exists(rdata_path.replace("raw.RData", coverage+".susie.RData")):
        os.system("rm %(rdata_path)s"%{"rdata_path":rdata_path})

###############################################

l = get_maf_list(outdir, coverage)
gwas_maf_set = l[0]
eqtl_maf_set = l[1]

gwas_rdata_set = set([i+".raw.RData" for i in gwas_maf_set])
eqtl_rdata_set = set([i+".raw.RData" for i in eqtl_maf_set])


pool = Pool(args.threads)
for rdata_path in gwas_rdata_set:
    gwas_id = rdata_path.split("/")[-2].replace(".GWAS_split","")
    gwas_type = gwas_dict[gwas_id]
    pool.apply_async(task_susie, (rdata_path, "gwas", gwas_type))
for rdata_path in eqtl_rdata_set:
    pool.apply_async(task_susie, (rdata_path, "eqtl", "NULL",))
pool.close()
pool.join()
