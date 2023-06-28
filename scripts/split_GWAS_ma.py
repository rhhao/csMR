# -*- coding: utf-8 -*-
#########################################################################
# File Name: split_GWAS_ma.py
# Author: Ruo-Han Hao and Feng Jiang
# Last Modified: 2023-01-08
# Description: 
# Usage:
#########################################################################
import os,sys,re
import argparse
parser = argparse.ArgumentParser(description="Split GWAS ma file into .maf files.")
parser.add_argument('-ma', '--ma_file', help="Input a ma file.", type=str, required=False)
parser.add_argument('-w', '--window_bp', help="window surrounding the index SNPs, default 100000bp.", type=int, default=100000, required=False)
parser.add_argument('-r', '--reference_genotype', help="Path to plink files of reference genotype.", type=str, required=True)
parser.add_argument('-id', '--gwas_id', help="GWAS ID.", type=str, required=True)
parser.add_argument('-o', '--outdir', help="Output directory.", type=str, required=True)

args = parser.parse_args()
ma_file = args.ma_file
if not os.path.isabs(ma_file):
    ma_file = os.path.join(os.getcwd(), ma_file)
window = args.window_bp
ref_path = args.reference_genotype
id = args.gwas_id
outdir = args.outdir

def trans_beta(snp, a1, a2, beta):
    bim = d.get(snp) #([(ch, pos, A1, A2),..])
#    bim_allele = d[snp][2:4]
#    if [a1, a2] == bim_allele:
    dic = {"A":"T", "T":"A", "G":"C", "C":"G"}
    for snp_info in bim:
        bim_allele = snp_info[2:4]
        if a1 == bim_allele[0] or [dic[i] if i in dic else i for i in [a1, a2]] == bim_allele:
            beta_new = beta
            break
#    elif [a2, a1] == bim_allele:
        elif a2 == bim_allele[0] or [dic[i] if i in dic else i for i in [a2, a1]] == bim_allele:
            beta_new =  "%.6g" % ((-1)*float(beta))
            break
        else:
            print(snp, [a1, a2], bim_allele, " allele not consistent with ref")
            beta_new = "NA"
    return beta_new


def task(ma_file, outdir, id):
    outdir = os.path.join(outdir, "precomputation", id+".GWAS_split")
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    os.system("ln -s %(ma_file)s %(outdir)s/%(ma_file_base)s" % {"ma_file": ma_file, "outdir": outdir, "ma_file_base": os.path.basename(ma_file)})
    allsnp_bedfile = outdir+"/"+id+".pos.bed"  # bed file include all SNPs
    sigsnp_bedfile = outdir+"/"+id+".5e-8SNP.100kb.bed"  # bed file include SNP with P<5e-8 and its 100kb flanking region
    d_ma = {}
    ma_file = ma_file
    with open(ma_file) as f1:
        firstline = f1.readline().strip().split("\t")
        indx_snp = firstline.index("SNP")
        indx_a1 = firstline.index("A1")
        indx_a2 = firstline.index("A2")
        indx_maf = firstline.index("MAF")
        indx_b = firstline.index("BETA")
        indx_s = firstline.index("SE")
        indx_p = firstline.index("P")
        indx_n = firstline.index("N")
        with open(allsnp_bedfile, 'w') as ff1:
            with open(sigsnp_bedfile, 'w') as ff2:
                for i in f1:
                    i = i.strip().split("\t")
                    i[indx_a1], i[indx_a2] = i[indx_a1].upper(), i[indx_a2].upper()
                    if 0 < float(i[indx_maf]) < 1 and float(i[indx_b]) != 0 and 0 < float(i[indx_p]) < 1:  # check if maf and p between 0,1 and beta != 0
                        i[indx_snp] = i[indx_snp].replace(":","_")
                        if i[indx_snp] in d and re.search("^[ACGT]+$", i[indx_a1]+i[indx_a2]):  # check if SNP was included in reference genotype and was coded by ACGT
                            i[indx_b] = trans_beta(i[indx_snp], i[indx_a1], i[indx_a2], i[indx_b]) # make sure the sign of beta is consistent with reference
                            if i[indx_b] != "NA":
                                i[indx_maf] = "%.6g" % (min(float(i[indx_maf]), 1-float(i[indx_maf])))
                                d_ma[i[indx_snp]] = [i[indx_snp], i[indx_maf], i[indx_b], i[indx_s], i[indx_p], i[indx_n]]  # snp:snp maf Beta    SE  P   N 
                                info = d.get(i[indx_snp])[0]
                                ch = info[0]
                                pos1 = str(int(info[1])-1)
                                pos2 = str(int(info[1])+1)
                                ff1.write("\t".join([ch, pos1, pos2, i[indx_snp]])+"\n")
                                if float(i[indx_p]) < 5e-8: 
                                    pos1 = str(max(int(info[1])-window, 0))
                                    pos2 = str(int(info[1])+window)
                                    ff2.write("\t".join([ch, pos1, pos2, i[indx_snp]])+"\n")
    merged_sigsnp_bedfile = outdir+"/"+id+".5e-8SNP.100kb.merged_nomhc.bed"
    os.system("bedtools subtract -a %(sigsnp_bedfile)s -b data/Genome_anno/MHC-extension-region-NCBI.nochr.bed |sort -k1,1 -k2,2n |bedtools merge -i - -c 4 -o collapse -delim \"_\" > %(merged_sigsnp_bedfile)s" %
              {"sigsnp_bedfile": sigsnp_bedfile, "merged_sigsnp_bedfile": merged_sigsnp_bedfile})
    d_maf = {}
    with os.popen('cut -f1-3 %(merged_sigsnp_bedfile)s|bedtools intersect -a - -b %(allsnp_bedfile)s -wo ' % {"merged_sigsnp_bedfile": merged_sigsnp_bedfile, "allsnp_bedfile": allsnp_bedfile}) as f1:
        for i in f1:
            i = i.strip().split("\t")
            d_maf.setdefault("_".join([i[0], i[1], i[2]]), set([])).add(i[6])  # ch_pos_pos:set([snp,snp,])
    for segment in d_maf:
        with open(outdir+"/"+segment+".maf", 'w') as ff:
            ff.write("SNP\tMAF\tBETA\tSE\tP\tN\n")
            for snp in d_maf[segment]:   
                ff.write("\t".join(d_ma[snp])+"\n")
#######################################################
d = {}
for ch in range(1, 23):
    pattern = "chr"+str(ch)+"\D+.*bim$"
    for file in os.listdir(ref_path):
        if re.search(pattern, file, re.I):
            bim_file = os.path.join(ref_path, file)
    with open(bim_file) as f1:
        for i in f1:
            i = i.strip().split("\t")
            i[1] = i[1].replace(":","_")
            d.setdefault(i[1],[]).append([i[0], i[3], i[4], i[5]])

task(ma_file, outdir, id)
for file in os.listdir(outdir+"/precomputation/"+id+".GWAS_split"):
    if re.search(".+\.pos.bed$", file) or re.search(".+\.5e-8SNP.100kb.bed$", file):
        os.remove(os.path.join(outdir, "precomputation", id+".GWAS_split", file))
