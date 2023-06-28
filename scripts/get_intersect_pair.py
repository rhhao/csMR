# -*- coding: utf-8 -*-
#########################################################################
# File Name: get_intersect_pair.py
# Author: Ruo-Han Hao and Feng Jiang
# Last Modified: 2023-01-08
# Description: find intersect pair between GWAS and eQTL
# Usage: python get_intersect_pair.py -gd [gwas_dir] -ed [eqtl_dir] -o [outdir]
#########################################################################
import os
import sys
import re
import argparse
import pybedtools
parser = argparse.ArgumentParser(description="Get intersection maf file pairs.")
parser.add_argument('-gd', '--gwas_split_dir', help="Output directory", type=str, required=True)
parser.add_argument('-ed', '--eqtl_split_dir', help="Output directory", type=str, required=True)
parser.add_argument('-o', '--outdir', help="Output directory", type=str, required=True)

args = parser.parse_args()
gwas_dir = args.gwas_split_dir
eqtl_dir = args.eqtl_split_dir
outdir = args.outdir

def get_intersect(eqtl_dir, gwas_dir, outdir):
    eqtl_maf_list = [i.split("--")[0].split("_")+[eqtl_dir+"/"+i] for i in os.listdir(eqtl_dir)
                     if re.search("\.maf$", i)]  # [[chr,pos1,pos2,eqtl_maf_file_path],]
    gwas_maf_list = [re.sub("\.maf$", "", i).split("_")+[gwas_dir+"/"+i]
                     for i in os.listdir(gwas_dir) if re.search("\.maf$", i)]  # [[chr,pos1,pos2,gwas_maf_file_path],]
    gwas_segment, eqtl_segment = pybedtools.BedTool(gwas_maf_list), pybedtools.BedTool(eqtl_maf_list)
    intersect_segment = gwas_segment.intersect(eqtl_segment, wo=True) 
    out_name = os.path.join(outdir, "precomputation", "coloc."+os.path.basename(gwas_dir)+"."+os.path.basename(eqtl_dir)+".pairs")
    with open(out_name, 'w') as ff:
        for i in intersect_segment:
            ff.write("\t".join([i[3], i[7], i[8]])+"\n")
    pybedtools.cleanup()

##################################################################
if gwas_dir != None and eqtl_dir != None:
    get_intersect(eqtl_dir, gwas_dir, outdir)
else:
    print("No GWAS directory or eQTL directory was found. Please check your output path.")
    sys.exit()
