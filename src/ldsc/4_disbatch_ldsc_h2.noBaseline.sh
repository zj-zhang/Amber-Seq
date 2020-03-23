#!/bin/bash
TARGET_DIR=$1
PARDIR="."
QUEUE=$2
if [ "$QUEUE" == "" ]; then
	QUEUE='ccb'
fi

fid="$(basename -- $TARGET_DIR)"
echo "$fid"
if [ -f 4_disBatch_$fid.noBaseline.txt ]; then
	rm 4_disBatch_$fid.noBaseline.txt
fi

for GWAS in `cat ./ldsc_resources/ukbb_sumstats_alkesprice/sumstats_list.txt`
do
    for ANNOT in `cat ./ldsc_resources/labels_annot/total_label_idx.csv | sed -e "1d" | awk '{print $(NF)}'`
    do
		    echo "./4_sbatch_ldsc_h5.noBaseline.sh $GWAS $ANNOT $PARDIR/$TARGET_DIR" >> 4_disBatch_$fid.noBaseline.txt
    done
done

sbatch -J 4_$fid -o slurm-4_ldsc_h5-%j.out -p $QUEUE --mem=500G --time=1-12 -N 2 -n 80 --ntasks-per-node 40 --wrap "disBatch.py 4_disBatch_$fid.noBaseline.txt"
