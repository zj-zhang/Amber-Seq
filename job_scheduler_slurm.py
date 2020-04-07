#!/usr/bin/env python3
"""
This script is originally from: https://groups.google.com/forum/#!topic/snakemake/7cyqAIfgeq4
I modified it for handling different GPU nodes.
ZZ
2020.3.30

Below is the original DOCSTRING:

Submit this clustering script for sbatch to snakemake with:

    snakemake -j 99 --cluster slurm_scheduler.py
"""

import os
import sys
import warnings


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

from snakemake.utils import read_job_properties

jobscript = sys.argv[1]

job_properties = read_job_properties(jobscript)

cluster_param={}

job_specs = job_properties["cluster"]


if not "mem" in job_specs:
    warnings.warn("Rule {rule} has no memory specified, set to default.".format(**job_properties))
    eprint(job_properties)

# do something useful with the threads
cluster_param['time'] = job_specs.get("time", job_properties["cluster"]['time'])
cluster_param['mem'] = job_specs.get("mem", '10GB') 
cluster_param['nodes'] = int(job_specs.get("num_nodes", 1) )
cluster_param['ntasks'] = int(job_specs.get("ntasks", 1) )
cluster_param['ntasks_per_node'] = int(job_specs.get("ntasks_per_node", 1))
cluster_param['queue'] = job_specs['queue']

cluster_param['name'] = job_properties['rule']

if 'gres' in job_specs:
    cluster_param['gres'] = job_specs['gres']
    use_gpu = True
    assert cluster_param['queue'] == 'gpu'
else:
    use_gpu = False
    assert cluster_param['queue'] != 'gpu'


# access property defined in the cluster configuration file (Snakemake >=3.6.0)
#job_properties["cluster"]["profile"]

eprint("Submit job with parameters:\n"+"\n".join(["\t{} : {}".format(key,cluster_param[key]) for key in cluster_param]))

if use_gpu:
    os.system("sbatch -p {queue} --gres={gres} -N {nodes} --ntasks {ntasks} --ntasks-per-node {ntasks_per_node} --time={time} --mem={mem} --job-name={name} {script}".format(script=jobscript, **cluster_param))
else:
    os.system("sbatch -p {queue} -N {nodes} --ntasks {ntasks} --ntasks-per-node {ntasks_per_node} --time={time} --mem={mem} --job-name={name} {script}".format(script=jobscript, **cluster_param))
