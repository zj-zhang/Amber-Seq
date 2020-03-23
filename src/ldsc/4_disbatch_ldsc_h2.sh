#!/bin/bash
TARGET_DIR=$1
PARDIR="."
QUEUE=$2
if [ "$QUEUE" == "" ]; then
	QUEUE='ccb'
fi

fid="$(basename -- $TARGET_DIR)"
echo "$fid"
if [ -f 4_disBatch_$fid.txt ]; then
	rm 4_disBatch_$fid.txt
fi

for GWAS in `cat ./ldsc_resources/ukbb_sumstats_alkesprice/sumstats_list.txt`
do
    for ANNOT in `cat ./ldsc_resources/labels_annot/total_label_idx.csv | sed -e "1d" | awk '{print $(NF)}'`
    do
		    echo "./4_sbatch_ldsc_h5.sh $GWAS $ANNOT $PARDIR/$TARGET_DIR" >> 4_disBatch_$fid.txt
    done
done

python -c "import BioNAS"
sbatch -J 4_$fid -o slurm-4_ldsc_h5-%j.out -p $QUEUE --mem=500G --time=3-12 -N 5 -n 200 --ntasks-per-node 40 --wrap "disBatch.py 4_disBatch_$fid.txt"
