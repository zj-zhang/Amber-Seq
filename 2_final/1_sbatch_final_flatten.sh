#!/bin/bash
#SBATCH -p gpu --gres=gpu:v100-32gb:1
##SBATCH --constraint="v100|p100"
#SBATCH --mem=64G -c 8
##SBATCH --exclusive
#SBATCH --time=3-13
##SBATCH -o slurm-finalConv-%x.out

# args: search-index layer train-index
sd=$1
od=$2
final=$4
wsf=$3
if [ ! -d $od ]; then
  mkdir -p $od
fi

FINAL_SCRIPT="convert_keras.py"
#FINAL_SCRIPT="final_conv.py"
TIMESTAMP=`date +"%m-%d-%y-%T"`


if [ "$final" == ".Random" ]; then
	python -u $FINAL_SCRIPT --gpus 1 --disable-controller --sd $sd --od $od --verbose 2 --wsf $wsf > $od/log.txt
else
	python -u $FINAL_SCRIPT --gpus 1 --sd $sd --od $od --verbose 2 --wsf $wsf > $od/log.txt
fi
