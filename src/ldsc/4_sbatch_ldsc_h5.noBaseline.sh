#!/bin/bash
#SBATCH -p ccb
#SBATCH --mem=16G -c 4
#SBATCH --time=2-12:00:00
##SBATCH -o s-%x.out

SUMSTATS=$1
ANNOT_LIST=$2
PAR_DIR=$3

sumstats=$(basename -- $SUMSTATS)
annot=$(basename -- $ANNOT_LIST)
if [ ! -d $PAR_DIR/vep_output/logs/ldsc_reg ]; then
	mkdir -p $PAR_DIR/vep_output/logs/ldsc_reg
fi	
python -u ldsc_reg_run_noBaseline.py $SUMSTATS $ANNOT_LIST $PAR_DIR 
