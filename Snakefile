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

configfile:
	CONFIG_FP


rule all:
	input:
		#"outputs/{project}/nas_final/metrics.log".format(project=PROJECT)
		"outputs/{project}/vep/".format(project=PROJECT) + config['vep']['vcf_prefix'] + "_abs_diffs.h5"

rule nas_search:
	input:
		CONFIG_FP
	output:
		"outputs/{project}/nas_search/train_history.png".format(project=PROJECT)

	params:
		controller_config_fp = config['nas_search']['controller_config_fp'],
		model_space_fp = config['nas_search']['model_space_fp'],
		working_dir = "outputs/{project}/nas_search/".format(project=PROJECT),
		pypath = PYTHONPATH
	log:
		"outputs/{project}/logs/search.log".format(project=PROJECT)

	shell:
		"export PYTHONPATH={params.pypath} \n"
		"python -u src/nas_search/search_conv.py --wd {params.working_dir} --model-space {params.model_space_fp} "
		"--controller-config {params.controller_config_fp} --verbose 2"


rule nas_final:
	input:
		"outputs/{project}/nas_search/train_history.png".format(project=PROJECT)

	output:
		"outputs/{project}/nas_final/metrics.log".format(project=PROJECT)

	params:
		final_config_fp = config['nas_final']['final_train_config_fp'],
		gpus = config['nas_final']['num_gpus'],
		width_scale_factor = config['nas_final']['width_scale_factor'],
		sd = "outputs/{project}/nas_search/".format(project=PROJECT),
		od = "outputs/{project}/nas_final/".format(project=PROJECT),
		pypath = PYTHONPATH

	log:
		"outputs/{project}/logs/final.log".format(project=PROJECT)

	shell:
		"export PYTHONPATH={params.pypath}\n"
		"python -u src/nas_final/convert_keras.py --config {params.final_config_fp} "
		"--gpus 1 "
	       	"--sd {params.sd} --od {params.od} --verbose 2 --wsf {params.width_scale_factor}"


rule variant_effect_prediction:
	input:
		"outputs/{project}/nas_final/metrics.log".format(project=PROJECT)

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
	input:
		"outputs/{project}/vep/" + config['vep']['vcf_prefix'] + "_abs_diffs.h5"
	
	output:
		"outputs/{project}/asb/{binding_type}/allelic_imbalance_summary.tsv"
	
	shell:
		"echo 'TODO'"
		

