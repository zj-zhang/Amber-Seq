#!/bin/bash
#SBATCH -p gpu --gres=gpu:1
#SBATCH --mem=128G -c 8
#SBATCH --time=3-13
#SBATCH -o s-vep-%x.out

# args: od
# first step: run vep
python -u vep.py --model $1/bestmodel.h5 --od $1/vep_output

# DEPRECATED: moved to 3_submit.sh
# second step: convert vep to ldsc annot
#python -u convert_ldsc.py $1/vep_output/ ./ldsc_resources/labels_annot/split_label_idx.$2.csv
