
########################################################################################
################################### FUNCTIONS ##########################################
########################################################################################
include: "rules/functions.smk"


########################################################################################
################################### VARIABLES ##########################################
########################################################################################
# Where all the output will be saved
import os,re
BASE_OUTPUT_DIR = config['BASE_OUTPUT_DIR']
if not os.path.exists(BASE_OUTPUT_DIR):
	os.makedirs(BASE_OUTPUT_DIR)

# Output file prefixes
RUN_PREFIXES = list(GWAS_SUMSTATS.keys())
CELL_TYPES = list(eQTL_SUMSTATS.keys())

# Defines reference genotype
GWAS_REFERENCE_PATH = config['GWAS_REFERENCE_GENOTYPE']['path']
if not os.path.isabs(GWAS_REFERENCE_PATH):
	GWAS_REFERENCE_PATH = os.path.join(os.getcwd(),GWAS_REFERENCE_PATH)
GWAS_REFERENCE_DUP = config['GWAS_REFERENCE_GENOTYPE']['duplicated_snp_path']
if not os.path.isabs(GWAS_REFERENCE_DUP):
	GWAS_REFERENCE_DUP = os.path.join(os.getcwd(),GWAS_REFERENCE_DUP)
eQTL_REFERENCE_PATH = config['eQTL_REFERENCE_GENOTYPE']
if not os.path.isabs(eQTL_REFERENCE_PATH):
	eQTL_REFERENCE_PATH = os.path.join(os.getcwd(),eQTL_REFERENCE_PATH)
	
# Set the max threads to do coloc analysis.
THREADS = config['COLOC_SETTING']['THREADS']

# coloc settings.
window_size = config['COLOC_SETTING']['WINDOW_SIZE_BP']
COVERAGES = str(config['COLOC_SETTING']['COVERAGES']).split(",")
COLOC_CUTOFF = config['COLOC_SETTING']['CUTOFF']

OUTCOME_DIR = config['OUTCOME_DIR']
OUTCOMRS = [os.path.splitext(file)[0] for file in os.listdir(OUTCOME_DIR)]

wildcard_constraints:
	coloc_cutoff = "0\.\d+",
	coverage = "0\.\d+"

########################################################################################
################################### Target files ##########################################
########################################################################################

list_target_files = []
#temp = expand("{BASE_OUTPUT_DIR}/{run_prefix}.GWAS_split/{run_prefix}.5e-8SNP.100kb.merged_nomhc.bed", BASE_OUTPUT_DIR = BASE_OUTPUT_DIR,
#run_prefix = RUN_PREFIXES)
#list_target_files.append(temp)
#temp = expand("{BASE_OUTPUT_DIR}/{run_prefix}.GWAS_split/{run_prefix}.GWAS_split.maflist",
#BASE_OUTPUT_DIR = BASE_OUTPUT_DIR, run_prefix = RUN_PREFIXES)
#list_target_files.append(temp)
#temp = expand("{BASE_OUTPUT_DIR}/{cell_type}.eQTL_split", BASE_OUTPUT_DIR = BASE_OUTPUT_DIR,
#cell_type = CELL_TYPES)
#list_target_files.append(temp)
#temp = expand("{BASE_OUTPUT_DIR}/coloc.{run_prefix}.GWAS_split.{cell_type}.eQTL_split.pairs", 
#BASE_OUTPUT_DIR = BASE_OUTPUT_DIR, run_prefix = RUN_PREFIXES, cell_type = CELL_TYPES)
#list_target_files.append(temp)
#temp = expand("{BASE_OUTPUT_DIR}/eqtl_maf.list", BASE_OUTPUT_DIR = BASE_OUTPUT_DIR)
#list_target_files.append(temp)
#temp = expand("{BASE_OUTPUT_DIR}/gwas_maf_finemap.list", BASE_OUTPUT_DIR = BASE_OUTPUT_DIR)
#list_target_files.append(temp)

#temp = expand("{BASE_OUTPUT_DIR}/eqtl_maf_finemap.list", BASE_OUTPUT_DIR = BASE_OUTPUT_DIR)
#list_target_files.append(temp)

#temp = expand("{BASE_OUTPUT_DIR}/coloc.{run_prefix}.GWAS_split.{cell_type}.eQTL_split.pairs.coloc_results", 
#BASE_OUTPUT_DIR = BASE_OUTPUT_DIR, run_prefix = RUN_PREFIXES, cell_type = CELL_TYPES)
#list_target_files.append(temp)

#temp = expand("{BASE_OUTPUT_DIR}/Susie_coloc.{run_prefix}.results.coverage.{coverage}",
#BASE_OUTPUT_DIR = BASE_OUTPUT_DIR, run_prefix = RUN_PREFIXES, coverage = COVERAGES)
#list_target_files.append(temp)

#temp = expand("{BASE_OUTPUT_DIR}/MR/Susie_coloc.{run_prefix}.results.coverage.{coverage}.PPH4.{coloc_cutoff}.{cell_type}.clumped",
#BASE_OUTPUT_DIR = BASE_OUTPUT_DIR, run_prefix = RUN_PREFIXES, coverage = COVERAGES, coloc_cutoff = COLOC_CUTOFF, cell_type = CELL_TYPES)
#list_target_files.append(temp)

temp = expand("{BASE_OUTPUT_DIR}/MR/{run_prefix}-{outcome}/{run_prefix}.{outcome}.{cell_type}.MR.res",
BASE_OUTPUT_DIR = BASE_OUTPUT_DIR, run_prefix = RUN_PREFIXES, cell_type = CELL_TYPES, outcome = OUTCOMRS)
list_target_files.append(temp)
temp = expand("{BASE_OUTPUT_DIR}/MR/{run_prefix}-{outcome}/{run_prefix}.{outcome}.{cell_type}.MR.sensitivity",
BASE_OUTPUT_DIR = BASE_OUTPUT_DIR, run_prefix = RUN_PREFIXES, cell_type = CELL_TYPES, outcome = OUTCOMRS)
list_target_files.append(temp)



########################################################################################
################################### PIPELINE ##########################################
########################################################################################
rule all:
	input:
		list_target_files

rule split_gwas:
	input:
		lambda wildcards: GWAS_SUMSTATS[wildcards.run_prefix]['path']
	output:
		"{BASE_OUTPUT_DIR}/precomputation/{run_prefix}.GWAS_split/{run_prefix}.5e-8SNP.100kb.merged_nomhc.bed"
	params:
		ma_file = lambda wildcards: GWAS_SUMSTATS[wildcards.run_prefix]['path'],
		id = "{run_prefix}",
		window_size = window_size,
		reference_genotype = GWAS_REFERENCE_PATH,
		outdir = BASE_OUTPUT_DIR
	#conda:
	#	"envs/envpy3.yml"
	shell:
		"python scripts/split_GWAS_ma.py -ma {params.ma_file} -w {params.window_size} -r {params.reference_genotype} -id {params.id} -o {params.outdir}"


rule gwas_rdata:
	input:
		"{BASE_OUTPUT_DIR}/precomputation/{run_prefix}.GWAS_split/{run_prefix}.5e-8SNP.100kb.merged_nomhc.bed"
	output:
		"{BASE_OUTPUT_DIR}/precomputation/{run_prefix}.GWAS_split/{run_prefix}.GWAS_split.maflist"
	params:
		reference_genotype = GWAS_REFERENCE_PATH,
		gwas_dir = "{BASE_OUTPUT_DIR}/precomputation/{run_prefix}.GWAS_split",
		gwas_type = lambda wildcards: GWAS_SUMSTATS[wildcards.run_prefix]['type'],
		threads = THREADS
#	conda:
#		"envs/envR4.yml"
	shell:
		"python scripts/extract_geno_ld.py -g {params.gwas_dir} -gr {params.reference_genotype} -gt {params.gwas_type} -t {params.threads}"


rule split_eqtl:
	input:
		expand("{BASE_OUTPUT_DIR}/precomputation/{run_prefix}.GWAS_split/{run_prefix}.GWAS_split.maflist",
			BASE_OUTPUT_DIR = BASE_OUTPUT_DIR, run_prefix = RUN_PREFIXES),
		lambda wildcards: eQTL_SUMSTATS[wildcards.cell_type]['path']
	output:
		directory("{BASE_OUTPUT_DIR}/precomputation/{cell_type}.eQTL_split")
	params:
		ma_file = lambda wildcards: eQTL_SUMSTATS[wildcards.cell_type]['path'],
		window_size = window_size,
		reference_genotype = eQTL_REFERENCE_PATH,
		outdir = BASE_OUTPUT_DIR
	#conda:
	#	"envs/envpy3.yml"
	shell:
		"python scripts/split_eQTL_ma.py -ma {params.ma_file} -w {params.window_size} -r {params.reference_genotype} -ct {wildcards.cell_type} -o {params.outdir}"


rule get_intersect_pair:
	input:
		"{BASE_OUTPUT_DIR}/precomputation/{cell_type}.eQTL_split"
	output:
		"{BASE_OUTPUT_DIR}/precomputation/coloc.{run_prefix}.GWAS_split.{cell_type}.eQTL_split.pairs"
	params:
		gwas_dir = "{BASE_OUTPUT_DIR}/precomputation/{run_prefix}.GWAS_split",
		eqtl_dir = "{BASE_OUTPUT_DIR}/precomputation/{cell_type}.eQTL_split",
		outdir = BASE_OUTPUT_DIR
	#conda:
	#	"envs/envpy3.yml"
	shell:
		"python scripts/get_intersect_pair.py -gd {params.gwas_dir} -ed {params.eqtl_dir} -o {params.outdir}"


rule eqtl_rdata:
	input:
		expand("{BASE_OUTPUT_DIR}/precomputation/coloc.{run_prefix}.GWAS_split.{cell_type}.eQTL_split.pairs",
				BASE_OUTPUT_DIR = BASE_OUTPUT_DIR, run_prefix = RUN_PREFIXES, cell_type = CELL_TYPES)
	output:
		"{BASE_OUTPUT_DIR}/precomputation/eqtl_maf.list"
	params:
		outdir = BASE_OUTPUT_DIR,
		reference_genotype = eQTL_REFERENCE_PATH,
		coverage = config['COLOC_SETTING']['COVERAGES']
#	conda:
#		"envs/envR4.yml"
	shell:
		"python scripts/extract_geno_ld.py --eqtl -o {params.outdir} -er {params.reference_genotype} -c {params.coverage}"


rule fine_mapping:
	input:
		"{BASE_OUTPUT_DIR}/precomputation/eqtl_maf.list"
	output:
		"{BASE_OUTPUT_DIR}/precomputation/gwas_maf_finemap.list",
		"{BASE_OUTPUT_DIR}/precomputation/eqtl_maf_finemap.list"
	params:
		outdir = BASE_OUTPUT_DIR,
		gwas_dict = GWAS_TYPE_DICT_STR,
		coverage = config['COLOC_SETTING']['COVERAGES'],
		threads = THREADS
#	conda:
#		"envs/envR4.yml"
	shell:
		"python scripts/susie_finemap.py -o {params.outdir} -d {params.gwas_dict} -c {params.coverage} -j {params.threads}"


rule run_coloc:
	input:
		expand("{BASE_OUTPUT_DIR}/precomputation/gwas_maf_finemap.list", BASE_OUTPUT_DIR = BASE_OUTPUT_DIR),
		expand("{BASE_OUTPUT_DIR}/precomputation/eqtl_maf_finemap.list", BASE_OUTPUT_DIR = BASE_OUTPUT_DIR)
	output:
		outdirs = directory([os.path.join(BASE_OUTPUT_DIR, "COLOC", "coloc.%s.GWAS_split.%s.eQTL_split.pairs.coloc_results"%(run_prefix, cell_type)) for run_prefix in RUN_PREFIXES for cell_type in CELL_TYPES ] )
	params:
		outdir = BASE_OUTPUT_DIR,
		coverage = config['COLOC_SETTING']['COVERAGES'],
		threads = THREADS
#	conda:
#		"envs/envR4.yml"
	shell:
		"python scripts/susie_coloc.py -o {params.outdir} -c {params.coverage} -t {params.threads}"

rule coloc_summary:
	input:
		expand("{{BASE_OUTPUT_DIR}}/COLOC/coloc.{{run_prefix}}.GWAS_split.{cell_type}.eQTL_split.pairs.coloc_results", cell_type = CELL_TYPES)
	output:
		"{BASE_OUTPUT_DIR}/COLOC/Susie_coloc.{run_prefix}.results.coverage.{coverage}.PPH4.{coloc_cutoff}"
	params:
		coverage = config['COLOC_SETTING']['COVERAGES'],
		id = "{run_prefix}",
		outdir = BASE_OUTPUT_DIR,
		coloc_cutoff = COLOC_CUTOFF
	#conda:
	#	"envs/envpy3.yml"
	shell:
		"python scripts/summarize.py -c {params.coverage} -id {params.id} -o {params.outdir} -co {params.coloc_cutoff}"

rule snp_filter:
	input:
		"{BASE_OUTPUT_DIR}/COLOC/Susie_coloc.{run_prefix}.results.coverage.{coverage}.PPH4.{coloc_cutoff}"
	output:
		"{BASE_OUTPUT_DIR}/MR/Susie_coloc.{run_prefix}.results.coverage.{coverage}.PPH4.{coloc_cutoff}.cells"
	params:
		coloc_out = "{BASE_OUTPUT_DIR}/COLOC/Susie_coloc.{run_prefix}.results.coverage.{coverage}.PPH4.{coloc_cutoff}",
		gwas_ma_file = lambda wildcards: GWAS_SUMSTATS[wildcards.run_prefix]['path'],
		cells = ",".join(CELL_TYPES),
		outdir = BASE_OUTPUT_DIR
#	conda:
#		"envs/envR4.yml"
	shell:
		"Rscript scripts/IV_select.r -r {params.coloc_out} -e {params.gwas_ma_file} -c {params.cells} -o {params.outdir}"

rule clump:
	input:
		"{BASE_OUTPUT_DIR}/MR/Susie_coloc.{run_prefix}.results.coverage.{coverage}.PPH4.{coloc_cutoff}.cells"
	output:
		"{BASE_OUTPUT_DIR}/MR/clump/Susie_coloc.{run_prefix}.results.coverage.{coverage}.PPH4.{coloc_cutoff}.{cell_type}.clumped"
	params:
		exp_ma_file = lambda wildcards: GWAS_SUMSTATS[wildcards.run_prefix]['path'],
		snp_file = "{BASE_OUTPUT_DIR}/MR/Susie_coloc.{run_prefix}.results.coverage.{coverage}.PPH4.{coloc_cutoff}.{cell_type}.IV",
		reference_genotype = GWAS_REFERENCE_PATH,
		reference_dupliaction_path = GWAS_REFERENCE_DUP
#	conda:
#		"envs/envR4.yml"
	shell:
		"Rscript scripts/clump.r -e {params.exp_ma_file} -s {params.snp_file} -r {params.reference_genotype} -d {params.reference_dupliaction_path}"

rule MR:
	input:
		expand("{{BASE_OUTPUT_DIR}}/MR/clump/Susie_coloc.{{run_prefix}}.results.coverage.{coverage}.PPH4.{coloc_cutoff}.{{cell_type}}.clumped",
		coverage = COVERAGES, coloc_cutoff = COLOC_CUTOFF)
	output:
		"{BASE_OUTPUT_DIR}/MR/{run_prefix}-{outcome}/{run_prefix}.{outcome}.{cell_type}.MR.res",
		"{BASE_OUTPUT_DIR}/MR/{run_prefix}-{outcome}/{run_prefix}.{outcome}.{cell_type}.MR.sensitivity"
	params:
		exp_ma_file = lambda wildcards: GWAS_SUMSTATS[wildcards.run_prefix]['path'],
		outcome_ma_file = OUTCOME_DIR+"/{outcome}.ma",
		clump_file = expand("{{BASE_OUTPUT_DIR}}/MR/clump/Susie_coloc.{{run_prefix}}.results.coverage.{coverage}.PPH4.{coloc_cutoff}.{{cell_type}}.clumped",
		coverage = COVERAGES, coloc_cutoff = COLOC_CUTOFF),
		ref_genotype = GWAS_REFERENCE_PATH
#	conda:
#		"envs/envR4.yml"
	shell:
		"Rscript scripts/csMR.r -e {params.exp_ma_file} -o {params.outcome_ma_file} -c {params.clump_file} -r {params.ref_genotype}"
