#!/bin/bash
PROJECT="$1" snakemake --cluster-config clusterconfig_slurm.yaml \
--cluster "sbatch -p {cluster.queue} -N {cluster.num_nodes} --ntasks {cluster.ntasks} --ntasks-per-node {cluster.ntasks_per_node} --mem={cluster.mem}" \
--latency-wait 60 --jobs 10 \
--rerun-incomplete \
$2
