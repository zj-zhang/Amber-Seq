#!/bin/bash
#GWAS="ldsc_resources/ukbb_sumstats/both_sexes/21001_irnt.ldsc.imputed_v3.both_sexes.tsv.gz"
#TARGET_DIR="tmp_baseline_DeepSEA/"
TARGET_DIR=$1
PARDIR="."

#PARDIR="./batch_run_20200117-L4_8_12-Flatten/"

#for TARGET_DIR in `ls $PARDIR | grep tmp_final`
#do
	for GWAS in `cat ./ldsc_resources/ukbb_sumstats_alkesprice/sumstats_list.txt`
	do
		for i in `seq 0 9`
		do
		    #sbatch -J $TARGET_DIR-$i ./4_sbatch_ldsc_h5.sh $GWAS ./ldsc_resources/labels_annot/split_label_idx.$i.csv $PARDIR/$TARGET_DIR
		    #bash 4_sbatch_ldsc_h5.sh $GWAS ./ldsc_resources/labels_annot/split_label_idx.$i.csv $PARDIR/$TARGET_DIR
		    echo "./4_sbatch_ldsc_h5.sh $GWAS ./ldsc_resources/labels_annot/split_label_idx.$i.csv $PARDIR/$TARGET_DIR"
		done
		#exit
	done
#done
