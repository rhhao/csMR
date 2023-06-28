# -*- coding: utf-8 -*-
#########################################################################
# File Name: susie_coloc.py
# Author: Ruo-Han Hao and Feng Jiang
# Last Modified: 2023-01-08
# Description: Calculate colocalization between eQTL and GWAS
# Usage: python susie_coloc.py -o [outdir] -c [coverage] -t [threads]
#########################################################################
import os
import sys
import re
from multiprocessing import Pool
import argparse
parser = argparse.ArgumentParser(description="Calculate colocalization between eQTL and GWAS")
#parser.add_argument('-p', '--pair_files', help="GWAS maf file path\teQTL maf file path, in a pair file. Multiple file split by comma", type=str, required=True)
parser.add_argument('-c', '--coverage', help="coverage in susie finemapping,  default 0.9, can split by comma", type=str, default=0.9, required=True)
parser.add_argument('-o', '--outdir', help="The output directory", type=str, required=True)
parser.add_argument('-t', '--threads', help="Threads number, default 3", type=int, default=3, required=False)

args = parser.parse_args()
#pair_files = args.pair_files
coverage_list = args.coverage.split(",")
outdir = args.outdir
threads = args.threads


def task(gwas_maf_path, eqtl_maf_path, out_dir):
    for coverage in coverage_list:
        if os.path.exists(gwas_maf_path+"."+coverage+".susie.RData") and os.path.exists(eqtl_maf_path+"."+coverage+".susie.RData"):
            os.system("Rscript scripts/susie_coloc_calculate.r -e %(gwas_maf_path)s -g %(eqtl_maf_path)s -o %(out_dir)s -c %(coverage)s" %
                    {"gwas_maf_path": gwas_maf_path, "eqtl_maf_path": eqtl_maf_path, "out_dir": out_dir, "coverage":coverage})
        else:
            if not os.path.exists(gwas_maf_path+"."+coverage+".susie.RData"):
                print(gwas_maf_path+"."+coverage+".susie.RData","not exists")
            if not os.path.exists(eqtl_maf_path+"."+coverage+".susie.RData"):
                print(eqtl_maf_path+"."+coverage+".susie.RData","not exists")
            continue

#################################
pool = Pool(threads)
coloc_pair = [file for file in os.listdir(os.path.join(outdir, "precomputation")) if re.search('^coloc\..*\.GWAS_split\..*\.eQTL_split\.pairs$', file)]
for pair_file in coloc_pair:
    out_dir = os.path.join(outdir, "COLOC", pair_file+".coloc_results")
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    pair_list = [i.strip().split("\t")[:2] for i in open(os.path.join(outdir, "precomputation", pair_file))]
    for pair in pair_list:
        pool.apply_async(task, (pair[0], pair[1], out_dir))
pool.close()
pool.join()
