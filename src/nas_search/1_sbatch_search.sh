#!/bin/bash
#SBATCH -p gpu --gres=gpu:v100-32gb:1
#SBATCH --mem=64G
#SBATCH --time=7-0
#SBATCH -o slurm-search-%x.out

parDir=$1
layer=$2
search=$3
if [ ! -d $parDir ]; then
	  mkdir -p $parDir
fi

DIR="$(cd "$(dirname "$0")" && pwd)"
python -u $DIR/search_conv.py --wd $parDir --layers $layer --verbose 2 > $parDir/log.txt
