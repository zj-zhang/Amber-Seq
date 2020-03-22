#!/bin/bash

parDir=$1
fid="$(basename -- $parDir)"
if [ -f $parDir/vep_output/hg19.baselineLD_All.20180423_abs_diffs.h5 ]; then
	echo "od=$parDir, fid=$fid"
	if [ -f 3_disBatch_$fid.txt ]; then
		rm 3_disBatch_$fid.txt
	fi
	for j in `cat ./ldsc_resources/labels_annot/split_chroms.txt`; do
		echo "bash 3_sbatch_ldsc_l2.sh $parDir $j" >> 3_disBatch_$fid.txt
	done
fi

cwd=`pwd`
python -c "import BioNAS"
sbatch -J 3_$fid -o slurm-%j.out -p ccb --mem=500G -c 5 --time=1-12 -N 1 -n 5 --ntasks-per-node 5 --exclusive --wrap "disBatch.py 3_disBatch_$fid.txt"
