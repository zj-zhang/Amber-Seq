#!/bin/bash
PROJECT="$1" snakemake --cluster-config clusterconfig_slurm.yaml \
--nolock \
--cluster ./job_scheduler_slurm.py \
--latency-wait 60 --jobs 10 \
--rerun-incomplete \
$2


#--cluster "sbatch -p {cluster.queue} -N {cluster.num_nodes} --ntasks {cluster.ntasks} --ntasks-per-node {cluster.ntasks_per_node} --mem={cluster.mem}" \
#--cluster-status ./job_status_slurm.py \
