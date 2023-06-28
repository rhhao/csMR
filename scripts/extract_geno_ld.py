#########################################################################
# File Name: extract_geno_ld.py
# Author: Ruo-Han Hao and Feng Jiang
# Last Modified: 2023-01-08
# Description: calculate LD matrix and convert maf, LD matrix into RData file.
# Usage:
#########################################################################
import os
import sys
import re
import argparse
from multiprocessing import Pool
parser = argparse.ArgumentParser(description="Extract genotpye and calculate ld matrix. Then generate RData. Note: when processing eQTL data, set sdY=1")
parser.add_argument('-g', '--gwas_maf_files_dir', help="Process all GWAS maf files in this directory", type=str, required=False)
parser.add_argument('-gr', '--gwas_reference_genotype', help="Path to plink files of GWAS reference genotype", type=str, required=False)
parser.add_argument('-gt', '--gwas_data_type', help="GWAS data type: cc or quant", type=str, required=False)
parser.add_argument('-e', '--eqtl', help="Process eQTL maf files", action = "store_true")
parser.add_argument('-o', '--outdir', help="The ouput directory", type=str, required=False)
parser.add_argument('-er', '--eqtl_reference_genotype', help="Path to plink files of eQTL reference genotype", type=str, required=False)
parser.add_argument('-c', '--coverage', help="Susie finemapping coverage", type=str, required=False)
parser.add_argument('-t', '--threads', help="Threads number, default 3", type=int, default=3, required=False)
args = parser.parse_args()
gwas_maf_dir = args.gwas_maf_files_dir
gwas_data_type = args.gwas_data_type
outdir = args.outdir
gwas_ref = args.gwas_reference_genotype
eqtl_ref = args.eqtl_reference_genotype
coverage = args.coverage
threads = args.threads

def get_eqtl_maf_list(outdir, coverage):
    outdir = os.path.join(outdir, "precomputation")
    pair_list = [os.path.join(outdir, file) for file in os.listdir(outdir) if re.search('^coloc\..*\.GWAS_split\..*\.eQTL_split\.pairs$', file)]
    with open(os.path.join(outdir, "eqtl_maf.list"), 'w') as fw:
        eqtl_maf_list=[]
        for i in pair_list:
            with open(i) as f:
                for line in f:
                    eqtl_maf = line.strip().split("\t")[1]
                    if not os.path.exists(eqtl_maf+"."+coverage+".susie.RData"):
                        eqtl_maf_list.append(eqtl_maf)
        eqtl_maf_list=set(eqtl_maf_list)
        fw.write("\n".join(eqtl_maf_list)+"\n")
    return(list(eqtl_maf_list))

def task(maf_file, ref_path, data_type, sdY_T):
    snp_file = re.sub("\.maf$", "", maf_file)+".snplist"
    with open(snp_file, 'w') as ff:
        with open(maf_file) as f1:
            next(f1)
            for i in f1:
                ff.write(i.strip().split("\t")[0]+"\n")
    ch = os.path.basename(maf_file).split("_")[0]
    pattern = "chr"+str(ch)+"\D+.*bim$"
    for file in os.listdir(ref_path):
        if re.search(pattern, file, re.I):
            geno_path = os.path.join(ref_path, file)
            geno_path = os.path.splitext(geno_path)[0]
    os.system("plink --silent --allow-no-sex --bfile %(geno_path)s --extract %(snp_file)s --make-bed --out %(snp_file)s" %
              {"geno_path": geno_path, "snp_file": snp_file})
    # ld calculation
    os.system("plink --silent --allow-no-sex --bfile %(snp_file)s --r bin4 --out %(snp_file)s" % {"snp_file": snp_file})
    os.system("rm %(snp_file)s %(snp_file)s.bed  %(snp_file)s.log %(snp_file)s.fam  %(snp_file)s.nosex" % {"snp_file": snp_file}) # keep bim file
    if sdY_T:
        os.system("Rscript scripts/process_raw2RData.r -m %(maf_file)s -t %(data_type)s -s TRUE" % {"maf_file": maf_file, "data_type": data_type})
    else:
        os.system("Rscript scripts/process_raw2RData.r -m %(maf_file)s -t %(data_type)s -s FALSE" % {"maf_file": maf_file, "data_type": data_type})
    ldbin_file = re.sub("\.maf$", "", maf_file)+".snplist.ld.bin"
    os.system("rm %(ldbin_file)s %(snp_file)s.bim" % {"ldbin_file": ldbin_file, "snp_file": snp_file})

#########################################################
if gwas_maf_dir != None:
    if gwas_ref == None:
        print("Please provide path to GWAS reference genotype!")
        sys.exit()
    else:
        os.system("ls %(gwas_maf_dir)s/*.maf > %(maf_file_list)s" % {"gwas_maf_dir": gwas_maf_dir, "maf_file_list": os.path.join(gwas_maf_dir,os.path.basename(gwas_maf_dir)+".maflist")})
        gwas_maf_set = set([gwas_maf_dir+"/"+i for i in os.listdir(gwas_maf_dir) if re.search("\.maf$", i)])
        pool = Pool(threads)
        for gwas_maf_file in gwas_maf_set:
            pool.apply_async(task, (gwas_maf_file, gwas_ref, gwas_data_type, False))
        pool.close()
        pool.join()
if args.eqtl:
    if outdir == None or eqtl_ref == None or coverage == None:
        print("Please check if you\n1)provide the output directory by -o/--outdir argument;\n2) provide the path to eQTL reference genotype by -er/--eqtl_reference_genotype argument;\n3)provide finemapping coverage you expected to use by -c/--coverage argument.")
        sys.exit()
    else:
        eqtl_maf_list = get_eqtl_maf_list(outdir, coverage)
        pool = Pool(threads)
        for eqtl_maf_file in eqtl_maf_list:
            pool.apply_async(task, (eqtl_maf_file, eqtl_ref, "quant", True))
        pool.close()
        pool.join()

