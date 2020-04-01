#!/bin/bash
PROJECT="$1" snakemake --cluster-config clusterconfig_slurm.yaml \
--cluster ./slurm_scheduler.py \
--latency-wait 60 --jobs 10 \
--rerun-incomplete \
$2


#--cluster "sbatch -p {cluster.queue} -N {cluster.num_nodes} --ntasks {cluster.ntasks} --ntasks-per-node {cluster.ntasks_per_node} --mem={cluster.mem}" \
#--cluster-status ./slurm_status.py \
