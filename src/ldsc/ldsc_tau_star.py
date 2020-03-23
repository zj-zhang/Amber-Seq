"""This script reads in necessary data and computes the standardized
effect sizes (i.e. tau^*) for LDscore regressions
ZZ
2020.2.25
"""

import os
import sys
import pandas as pd
import numpy as np
import h5py
import pickle
from tqdm import tqdm
from utils import run_from_ipython


def read_log(log_fp):
    #log_fp = os.path.join(wd, "ldsc_reg.log")
    total_snps = None
    h2 = None
    mean_chi2 = None
    lambda_gc = None
    intercept = None
    ratio = None
    with open(log_fp, "r") as f:
        for line in f:
            ele = line.strip().split()
            if line.startswith("Removed"):
                total_snps = int(ele[-3].strip("("))
            elif line.startswith("Total Observed scale h2"):
                h2 = float(ele[-2])
            elif line.startswith("Mean Chi"):
                mean_chi2 = float(ele[-1])
            elif line.startswith("Ratio"):
                ratio = float(ele[-2])
            elif line.startswith("Lambda GC"):
                lambda_gc = float(ele[-1])
            elif line.startswith("Intercept"):
                intercept = float(ele[-2])
            else:
                continue
    metrics ={"total_snps": total_snps, 
            "h2": h2, 
            "mean_chi2": mean_chi2, 
            "lambda_gc": lambda_gc, 
            "intercept": intercept, 
            "ratio": ratio
            }
    return metrics


def read_reg_coef(res_fp):
    #res_fp = os.path.join(wd, "ldsc_reg.results")
    df = pd.read_csv(res_fp, sep="\t", header=0)
    tau = df.iloc[0]['Coefficient']
    return tau


def get_annot_sd(pred_fp):
    par_dir = os.path.dirname(pred_fp)
    precompute_file = os.path.join(par_dir, "annot_std_list.pkl") 
    if os.path.isfile(precompute_file ):
        print("load precompute")
        sds = pickle.load(open(precompute_file, "rb"))
    else:
        if pred_fp.endswith("h5"):
            print("load h5")
            store = h5py.File(pred_fp, "r")
            pred_all = store['data'].value
            store.close()
        elif pred_fp.endswith("npz"):
            print("load npz")
            pred_all = np.load(pred_fp)['arr_0']
        else:
            raise Exception("Unknown data format: %s"%pred_fp)
        sds = np.apply_along_axis(np.std, axis=0, arr=pred_all)
        del pred_all
        pickle.dump(sds, open(precompute_file, "wb"))
    return sds


def get_tau_for_model(model_wd):
    vep_dir = "vep_output"
    #vep_dir = "vep_output_logfc" # this used to be "vep_output"
    score_fn = "hg19.baselineLD_All.20180423_abs_diffs.h5"
    #score_fn = "VEP_abs_logfc.npz"  # this used to be "hg19.baselineLD_All.20180423_abs_diffs.h5"
    sds = get_annot_sd(os.path.join(model_wd, vep_dir, score_fn))
    
    reg_dir = "label_wise_ldsc_reg_noBaseline_withLDblock" # this used to be "label_wise_ldsc_reg"  
    gwas_list = os.listdir( os.path.join(model_wd, vep_dir, reg_dir) )
    for gwas in gwas_list:
        print(gwas)
        annot_list = os.listdir( os.path.join(model_wd, vep_dir, reg_dir, gwas) )
        #total_snps, h2 = read_log( os.path.join(model_wd, vep_dir, reg_dir,gwas,  annot_list[0], "ldsc_reg.log") )
        try:
            for annot in tqdm(annot_list):
                tau_c = read_reg_coef( os.path.join(model_wd, vep_dir, reg_dir, gwas, annot, "ldsc_reg.results") )
                metrics = read_log( os.path.join(model_wd, vep_dir, reg_dir, gwas, annot, "ldsc_reg.log") )
                total_snps, h2 = metrics['total_snps'], metrics['h2']
                index = int(annot.split('--')[0])
                std = sds[index]
                tau_star = tau_c * total_snps / h2 * std

                output_fp =  os.path.join(model_wd, vep_dir, reg_dir, gwas, annot, "tau_star.txt")
                with open(output_fp, "w" ) as fo:
                    fo.write("\t".join(["tau_star", "tau_c", "total_snps", "h2", "std"]) + "\n")
                    fo.write("\t".join([str(x) for x in [tau_star, tau_c, total_snps, h2, std]] ) + "\n" )
        except Exception as e:
            print(e)
            break

if __name__ == '__main__':
    if not run_from_ipython():
        wd = sys.argv[1]
        if wd:
            get_tau_for_model(wd)




