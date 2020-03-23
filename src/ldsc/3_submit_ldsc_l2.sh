#!/bin/bash
#parDir="./batch_run_20200117-L4_8_12-Flatten"

#fileList=`ls $parDir | grep "tmp_final_*.*" | grep -v "tmp_final_12_5.1"`

#for parDir in $fileList
#do
parDir=$1
  if [ -f $parDir/vep_output/hg19.baselineLD_All.20180423_abs_diffs.h5 ]; then
    echo "od=$parDir"
    for j in `cat ./ldsc_resources/labels_annot/split_chroms.txt`; do
      sbatch -J ldsc-$i-$j 3_sbatch_ldsc_l2.sh $parDir $j
      #echo "bash 3_sbatch_ldsc_l2.sh $parDir/$i $j"
      
    done
  fi
#done
