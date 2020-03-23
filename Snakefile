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
		"outputs/{project}/nas_final/loss.png".format(project=PROJECT)


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
		"export PYTHONPATH={params.pypath}"
		"python -u src/nas_search/search_conv.py --wd {params.working_dir} --model-space {params.model_space_fp} "
		"--controller-config {params.controller_config_fp} --verbose 2"


rule nas_final:
	input:
		"outputs/{project}/nas_search/train_history.png".format(project=PROJECT)

	output:
		"outputs/{project}/nas_final/loss.png".format(project=PROJECT)

	params:
		gpus = config['nas_final']['num_gpus'],
		width_scale_factor = config['nas_final']['width_scale_factor'],
		sd = "outputs/{project}/nas_search/".format(project=PROJECT),
		od = "outputs/{project}/nas_final/".format(project=PROJECT),
		pypath = PYTHONPATH

	log:
		"outputs/{project}/logs/final.log".format(project=PROJECT)

	shell:
		"export PYTHONPATH={params.pypath}\n"
		"python -u src/nas_final/convert_keras.py --gpus 1 "
	       	"--sd {params.sd} --od {params.od} --verbose 2 --wsf {params.width_scale_factor}"


