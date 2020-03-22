#!/bin/bash
#parDir="./batch_run_20200117-L4_8_12-Flatten"
parDir=$1
fileList=`ls $parDir | grep tmp_final_12_[12] | grep -v Random`
for i in $fileList
do
  if [ -f $parDir/$i/loss.png ]; then
    echo "od=$i"
    echo "sbatch -J vep-$i-$j sbatch_vep.sh $parDir/$i"
    sbatch -J vep-$i-$j 2_sbatch_vep.sh $parDir/$i
  fi
done
