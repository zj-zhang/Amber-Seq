#!/bin/bash
parDir=$1
mode=$2
wsf=$3
if [ "$mode" == "nosearch" ]; then
	search=".noSearch"
	final=".noSearch"
elif [ "$mode" == "random" ]; then
	search=""
	final=".Random"
else
	search=""
	final=""
fi

if [ "$wsf" == "" ]; then
	wsf="2"
fi

for index in `seq 1 6`
do
  for layer in `seq 12 4 12`
  do
	sd=$parDir/tmp_search_"$layer"_"$index""$search"
	if [ "$wsf" == "2" ]; then
		od=$parDir/tmp_final_"$layer"_"$index""$final"
	else
		od=$parDir/tmp_final_"$layer"_"$index""$final"_wsf"$wsf"
	fi
	echo "sd=$sd, od=$od"	
	sbatch -J cf$index-$layer-$final 1_sbatch_final_flatten.sh $sd $od $wsf $final
  done
done
