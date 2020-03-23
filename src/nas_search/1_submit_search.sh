#!/bin/bash
parDir=$1

DIR="$(cd "$(dirname "$0")" && pwd)"
if [ ! -d $parDir ]
then
	mkdir -p $parDir
fi

for index in `seq 1 6`
do
  for layer in `seq 12 4 12`
  do
    echo "index=$index, layer=$layer"
    sbatch -J cs$index-$layer $DIR/1_sbatch_search.sh $parDir/bionas_search_"$layer"_"$index" $layer
  done
done
