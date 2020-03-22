#!/bin/bash
#SBATCH -p ccb
#SBATCH --mem=100G -c 1
#SBATCH --time=1-12:00:00
#SBATCH -o slurm-ldsc_l2-%j.out

# args: od chroms
OD=$1
CHROMS=$2
echo $OD
echo $CHROMS
VEP_DIR="vep_output" # used to be "vep_output"

if [ ! -d $OD/$VEP_DIR/logs/ ]; then
	mkdir -p $OD/$VEP_DIR/logs/
fi	
python -u convert_ldsc.py $OD/$VEP_DIR/ $CHROMS > $OD/$VEP_DIR/logs/$CHROMS.ldsc_annot.log 2>&1
#python -u convert_ldsc.py $OD/$VEP_DIR/ $CHROMS


