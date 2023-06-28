# -*- coding: utf-8 -*-
#########################################################################
# File Name: split_eQTL_ma.py
# Author: Ruo-Han Hao and Feng Jiang
# Last Modified: 2023-01-08
# Description: 
# Usage: 
#########################################################################
import os, sys, re, gc
from multiprocessing import Pool
import argparse
import pybedtools
parser = argparse.ArgumentParser(description="Split eQTL ma file into .maf files. File format is same as GTEX.")
parser.add_argument('-ma', '--ma_file', help="Input a ma file.", type=str, required=False)
parser.add_argument('-w', '--window_bp', help="window surrounding the index SNPs, default 100000bp.", type=int, default=100000, required=False)
parser.add_argument('-r', '--reference_genotype', help="Path to plink files of reference genotype.", type=str, required=True)
parser.add_argument('-ct', '--cell_type', help="eQTL cell type.", type=str, required=False)
parser.add_argument('-o', '--out_dir', help="Output directory.", type=str, required=False)
parser.add_argument('-t', '--threads', help="Threads used, default = 1.", type=int, default = 1)
args = parser.parse_args()
ma_file = args.ma_file
if not os.path.isabs(ma_file):
    ma_file = os.path.join(os.getcwd(), ma_file)
window=args.window_bp
ref_path = args.reference_genotype
cell_type = args.cell_type
outdir =  args.out_dir
threads = args.threads

def trans_beta(snp, a1, a2, beta):
    bim = d_snp.get(snp) #([(ch, pos, A1, A2),..])
    #bim_allele = d_snp[snp][2:4]
    dic = {"A":"T", "T":"A", "G":"C", "C":"G"}
    for snp_info in bim:
        bim_allele = snp_info[2:4]
        if a1 == bim_allele[0] or [dic[i] if i in dic else i for i in [a1, a2]] == bim_allele:
#    if [a1, a2] == bim_allele:
            beta_new = beta
            break    
        elif a2 == bim_allele[0] or [dic[i] if i in dic else i for i in [a2, a1]] == bim_allele:
            beta_new =  "%.6g" % ((-1)*float(beta))
            break
#    elif [a2, a1] == bim_allele:
        else:
            print(snp, [a1, a2], bim_allele, " allele not consistent with ref")
            beta_new = "NA"
    return beta_new

def get_d(ma_file): 
    d_ma = {}
    d_pos={}
    d_p={}
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
        indx_gene = firstline.index("GENE")
        for i in f1:
            i = i.strip().split("\t")
            if 0 < float(i[indx_maf]) < 1 and float(i[indx_b]) != 0 and 0 < float(i[indx_p]) < 1:  # check if maf and p between 0,1 and beta != 0
                i[indx_snp] = i[indx_snp].replace(":","_") 
                if re.search("^[ACGT]+$", i[indx_a1]+i[indx_a2]) and i[indx_snp] in d_snp:  # check if SNP was included in reference genotype and was coded by ACGT
                    i[indx_b] = trans_beta(i[indx_snp], i[indx_a1], i[indx_a2], i[indx_b])  # make sure the sign of beta is consistent with reference
                    if i[indx_b] != "NA":
                        i[indx_maf] = "%.6g" % (min(float(i[indx_maf]), 1-float(i[indx_maf])))
                        d_ma.setdefault(i[indx_gene], {})[i[indx_snp]]=[i[indx_snp], i[indx_b], i[indx_s], i[indx_p], i[indx_maf], i[indx_n]]  # ensg: {SNP:[SNP   b   se  p   MAF    N]}
                        snp_chr,snp_pos=d_snp[i[indx_snp]][0][:2]
                        snp_pos=int(snp_pos)
                        d_pos.setdefault(i[indx_gene],{})[i[indx_snp]]=[snp_chr,snp_pos,snp_pos,i[indx_snp]] #ensg:{rsid:[ chr pos pos rsid],}
                        p=float(i[indx_p])
                        if p<=0.05:
                            d_p.setdefault(i[indx_gene],{}).setdefault(p,set([])).add(i[indx_snp])  #ensg:{p:set([rs,rs,]),}

    d_top={} 
    for ensg in d_p:
        p_min=min(d_p[ensg])
        topsnp_set=d_p[ensg][p_min]
        d_top[ensg]=topsnp_set #ensg:set([rs,rs])

    return d_ma,d_pos,d_top

def extract_ensg_signal(ensg):
    if ensg in d_top: 
        sig_bed=[d_pos[ensg][snp] for snp in d_top[ensg]]
    else:
        print("No snp was found associated (P < 0.05) with", ensg)
        return
    sig_bed=[[i[0],max(i[1]-window,0),i[2]+window,i[3]] for i in sig_bed] 
    sig_bed=sorted(sig_bed,key=lambda x:x[1]) 
    sig_bed=pybedtools.BedTool(sig_bed)
    merge_sig_bed=sig_bed.merge()
    for segment in merge_sig_bed:
        ma_snplist=[snp for snp in d_ma[ensg]]
        ma_bed=pybedtools.BedTool([d_pos[ensg][snp] for snp in ma_snplist])
        segment_bed=pybedtools.BedTool([list(segment)])
        intersect_bed=ma_bed.intersect(segment_bed)
        out_name="_".join(segment)+"--"+ensg+"--"+cell_type+".maf"
        with open(out_dir+"/"+out_name,'w') as ff:
            ff.write("SNP\tBETA\tSE\tP\tMAF\tN\n")
            for i in intersect_bed:
                ff.write("\t".join(d_ma[ensg][i[3]])+"\n")
        pybedtools.cleanup()
    pybedtools.cleanup()

###################################################
d_snp = {}
for ch in range(1, 23): 
    pattern="chr"+str(ch)+"\D+.*bim$"
    for file in os.listdir(ref_path):
        if re.search(pattern, file, re.I):
            bim_file = os.path.join(ref_path, file)
    with open(bim_file) as f1:
        for i in f1:
            i = i.strip().split("\t")
            i[1] = i[1].replace(":","_")
            d_snp.setdefault(i[1],[]).append([i[0], i[3], i[4], i[5]])# rs:[chr,pos,A1,A2

out_dir=os.path.join(outdir, "precomputation", cell_type+".eQTL_split")
if not os.path.exists(out_dir):
    os.mkdir(out_dir)
with open(ma_file) as f:
    headline = f.readline()
    gene_indx = headline.strip().split("\t").index("GENE")
    genel = set([])
    sub_ma = ma_file+".sub"
    fw = open(sub_ma,"w")
    fw.write(headline)
    for line in f:
        i = line.strip().split("\t")
        genel.add(i[gene_indx])
        if len(genel) <= 2000:
            fw.write(line)
        else:
            fw.close()
            d_ma,d_pos,d_top=get_d(sub_ma)
            pool=Pool(threads)
            for ensg in d_pos:
                pool.apply_async(extract_ensg_signal,(ensg,))
            pool.close()
            pool.join()
            fw = open(sub_ma,"w")
            fw.write(headline)
            fw.write(line)
            genel = set([])
            genel.add(i[gene_indx])
    if len(genel) != 0:
        fw.close()
        d_ma,d_pos,d_top=get_d(sub_ma)
        d_snp = None
        gc.collect()
        pool=Pool(threads)
        for ensg in d_pos:
            pool.apply_async(extract_ensg_signal,(ensg,))
        pool.close()
        pool.join()
os.remove(sub_ma)

