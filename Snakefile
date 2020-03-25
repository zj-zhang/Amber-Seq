"""
A pipeline for optimizing DeepSEA architectures using BioNAS
Each `$PROJECT` is configured by its own `$CONFIG_FP`

Author: 
	Zijun Zhang <zj.z@ucla.edu>

Date: 
	3.22.2020

Histroy:

"""


import os
import sys


# IN-FILE CONFIGURATION
PROJECT = os.environ.get("PROJECT")
assert PROJECT is not None and len(PROJECT)>0, "must provide `$PROJECT`, e.g. `PROJECT=test snakemake -n`"
CONFIG_FP = os.path.join("outputs", PROJECT, 'config', 'config.yaml')
PYTHONPATH = os.path.abspath(os.getcwd())
print("python path: ", PYTHONPATH)
CHROM_PARTS = ["1,10,11,20", "2,9,12,19", "3,8,13,18", "4,7,14,17,22", "5,6,15,16,21"]

configfile:
	CONFIG_FP


rule all:
	input:
		#"outputs/{project}/vep/".format(project=PROJECT) + config['vep']['vcf_prefix'] + "_abs_diffs.h5"
		["outputs/{project}/asb/{binding_type}/allelic_imbalance_summary.tsv".format(project=PROJECT, binding_type=x) for x in config['allelic_imbalance']],
		"outputs/{project}/ldsc/label_wise_l2/done.txt".format(project=PROJECT),
		"outputs/{project}/ldsc/label_wise_h2/done.txt".format(project=PROJECT) 




rule nas_search:
	input:
		"outputs/{project}/config/config.yaml"
	output:
		"outputs/{project}/nas_search/train_history.png"

	params:
		controller_config_fp = config['nas_search']['controller_config_fp'],
		model_space_fp = config['nas_search']['model_space_fp'],
		working_dir = "outputs/{project}/nas_search/",
		pypath = PYTHONPATH
	log:
		"outputs/{project}/logs/search.log"

	shell:
		"export PYTHONPATH={params.pypath} \n"
		"python -u src/nas_search/search_conv.py --wd {params.working_dir} --model-space {params.model_space_fp} "
		"--controller-config {params.controller_config_fp} --verbose 2"


rule nas_final:
	input:
		"outputs/{project}/nas_search/train_history.png"

	output:
		"outputs/{project}/nas_final/metrics.log"

	params:
		final_config_fp = config['nas_final']['final_train_config_fp'],
		gpus = config['nas_final']['num_gpus'],
		width_scale_factor = config['nas_final']['width_scale_factor'],
		sd = "outputs/{project}/nas_search/",
		od = "outputs/{project}/nas_final/",
		pypath = PYTHONPATH

	log:
		"outputs/{project}/logs/final.log"

	shell:
		"export PYTHONPATH={params.pypath}\n"
		"python -u src/nas_final/convert_keras.py --config {params.final_config_fp} "
		"--gpus 1 "
	       	"--sd {params.sd} --od {params.od} --verbose 2 --wsf {params.width_scale_factor}"


rule variant_effect_prediction:
	input:
		"outputs/{project}/nas_final/metrics.log"

	output:
		"outputs/{project}/vep/" + config['vep']['vcf_prefix'] + "_abs_diffs.h5"

	params:
		model_path = "outputs/{project}/nas_final/bestmodel.h5",
		output_dir = "outputs/{project}/vep/",
		genome_fp = config['vep']['genome_fp'],
		vcf_fp = config['vep']['vcf_fp'],
		label_annot_fp = config['vep']['label_annot_fp'],
		pypath = PYTHONPATH
	
	shell:
		"export PYTHONPATH={params.pypath}\n"
		"python src/asb/vep.py --model {params.model_path} "
		"--od {params.output_dir} "
		"--genome {params.genome_fp} "
		"--vcf {params.vcf_fp} "
		"--label-annot {params.label_annot_fp} "


rule allelic_imbalance_analysis:
	# NOTE: for now, this rule relies on a pre-existing file that maps each 
	# Allelic imbalance SNP to a pair of Col/Row index in the h5py file
	input:
		"outputs/{project}/vep/" + config['vep']['vcf_prefix'] + "_abs_diffs.h5"
	
	output:
		"outputs/{project}/asb/{binding_type}/allelic_imbalance_summary.tsv"
	
	params:
		vep_prefix = "outputs/{project}/vep/" + config['vep']['vcf_prefix'],
		ref_pkl = lambda wildcards: config["allelic_imbalance"][wildcards.binding_type]['ref_pkl']

	shell:
		"export PYTHONPATH="+PYTHONPATH+"\n"
		"python src/asb/allelic_pred.py --vep-prefix {params.vep_prefix} "
		"--ref-pkl {params.ref_pkl} "
		"--out  {output} "


rule ldsc_l2_step:
	# NOTE: in order to use disbatch, the command-line are written to a text file;
	# to just run it, remove the echo at the header
	input:
		"outputs/{project}/vep/" + config['vep']['vcf_prefix'] + "_abs_diffs.h5"
	
	output:
		"outputs/{project}/ldsc/label_wise_l2/done.txt"
	
	params:
		cmd_list = "outputs/{project}/ldsc/label_wise_l2/ldsc_l2.cmd_list.txt",
		chrom_parts_fp = config["ldsc"]["chrom_parts_fp"],
		vep_dir = "outputs/{project}/vep/",
		label = config['ldsc']['label_fp'],
		ldsc_bfile = config['ldsc']['bfile'],
		ldsc_snps = config['ldsc']['snps'],
		ldsc_baselineLD = config['ldsc']['baselineLD'],
		ldsc_bin = config['ldsc']['bin'],
		output_prefix = "outputs/{project}/ldsc/label_wise_l2/disbatch_run_l2"	
	
	shell:
		"if [ -f {params.cmd_list} ]; then rm {params.cmd_list}; fi \n"
		"for chrom in `cat {params.chrom_parts_fp}`;"
		"do \n"
		"echo python src/ldsc/convert_ldsc.py --vep-dir {params.vep_dir} "
		"--chroms $chrom "
		"--label {params.label} "
		"--ldsc-bfile {params.ldsc_bfile} "
		"--ldsc-snps {params.ldsc_snps} "
		"--ldsc-baselineLD {params.ldsc_baselineLD} "
		"--ldsc-bin {params.ldsc_bin} "
		" >> {params.cmd_list} \n"
		"done"
		"disBatch.py -p {params.output_prefix} {params.cmd_list}\n"
		"echo `date` Done > {output}"


rule ldsc_h2_step:
	input:
		"outputs/{project}/ldsc/label_wise_l2/done.txt"
	
	output:
		"outputs/{project}ldsc/label_wise_h2/done.txt"

	params:
		cmd_list = "outputs/{project}/ldsc/label_wise_h2/ldsc_h2.cmd_list.txt",
		par_dir = "outputs/{project}/",
		sumstats_fp = config['ldsc']['sumstats'],
		annot_fp = config['ldsc']['label_fp'],
		ldsc_baselineLD = config['ldsc']['baselineLD'],
		ldsc_frqfile = config['ldsc']['frqfile'],
		ldsc_weight = config['ldsc']['weight'],
		ldsc_bin = config['ldsc']['bin'],
		ldsc_quantile_M = config['ldsc']['quantile_M'],
		ldsc_quantile_h2g = config['ldsc']['quantile_h2g'],
		output_prefix = "outputs/{project}/ldsc/label_wise_h2/disbatch_run_h2"

	shell:
		"if [ -f {params.cmd_list} ]; then rm {params.cmd_list}; fi \n"
		"for GWAS in `cat {params.sumstats_fp}`; do \n"
		"for ANNOT in `cat {params.annot_fp}`; do \n"
		"echo python src/ldsc/ldsc_reg_run.py --sumstats $GWAS"
		"--annot-list $annot "
		"--par-dir {params.par_dir} "
		"--ldsc-frqfile {params.ldsc_frqfile} "
		"--weight {params.ldsc_weight} "
		"--ldsc-quantile-M {params.ldsc_quantile_M} "
		"--ldsc-quantile-h2g {params.ldsc_quantile_h2g} "
		"--ldsc-bin {params.ldsc_bin} "
		">> {params.cmd_list} \n"
		"done\n done\n "
		"disBatch.py -p {params.output_prefix} {params.cmd_list} \n"
		"echo `date` Done > {output}"

