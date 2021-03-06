import subprocess
import sys
import os
from src.utils import deepsea_label_name_normalizer

import argparse

parser = argparse.ArgumentParser(description="Step 2 of LDSC: run regression on l2")
parser.add_argument("--sumstats", type=str, help="filepath to GWAS summary statistics")
parser.add_argument("--annot-list", type=str, help="A comma-separated string of Annotations to work on")
parser.add_argument("--par-dir", type=str, help="filepath to parent directory (i.e. project folder)")
parser.add_argument("--ldsc-frqfile", type=str, help="filepath to LDSC frequency file")
parser.add_argument("--baselineLD", type=str, help="filepath to baselineLD")
parser.add_argument("--weight", type=str, help="filepath to LDSC SNP weights")
parser.add_argument("--ldsc-quantile-M", type=str, help="filepath to executable quantile_M.pl")
parser.add_argument("--ldsc-quantile-h2g", type=str, help="filepath to executable quantile_h2g.r")
parser.add_argument("--ldsc-bin", type=str, help="filepath to executable (chmod +x) LDSC python script")

args = parser.parse_args()

sumstats = args.sumstats
annot_list = args.annot_list
par_dir = args.par_dir

resume_prev_run = True
#ldsc_bin = "./ldsc_resources/ldsc/ldsc.py"
ldsc_bin = args.ldsc_bin
#ldsc_frqfile =  "./ldsc_resources/frq_file/1000G_Phase3_frq/1000G.EUR.QC."
ldsc_frqfile = args.ldsc_frqfile
#ldsc_baselineLD = "./ldsc_resources/baselineLD/baselineLD."
ldsc_baselineLD = args.baselineLD
#ldsc_weight =  "./ldsc_resources/weight_file/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC."
ldsc_weight = args.weight
#ldsc_quantile_M_bin = "./ldsc_resources/ldsc/ContinuousAnnotations/quantile_M.pl"
ldsc_quantile_M_bin = args.ldsc_quantile_M
#ldsc_quantile_h2g_bin = "./ldsc_resources/ldsc/ContinuousAnnotations/quantile_h2g.r"
ldsc_quantile_h2g_bin = args.ldsc_quantile_h2g

sumstats_prefix = os.path.basename(sumstats).split(".")[0]
vep_dir = "label_wise_l2" # this used to be "vep_output"
reg_subdir = "label_wise_h2" #  "label_wise_ldsc_reg_BaselineV2.2"
reg_dir = os.path.join(par_dir, reg_subdir, sumstats_prefix)
os.makedirs(reg_dir, exist_ok=True)


if os.path.isfile(annot_list):
    # DEPRECATED: parsed a filename for a set of annots
    annot_list = [x.split("\t")[0] for x in open(annot_list,'r').readlines()]
    annot_list.pop(0)
else:
    # NEW: direct parse comma-separated annots
    annot_list = annot_list.split(",")

# 1. run ldsc-reg for each annot
print("1. run ldsc-reg for each annot")
for annot in annot_list:
    annot = deepsea_label_name_normalizer(annot)
    annot_fp = os.path.join(par_dir, vep_dir, "label_wise_ldsc", "%sL2"%annot, "%sL2."%annot)
    if not os.path.isfile(annot_fp+"1.l2.ldscore.gz"):
        print("%s File not found: %s"%(annot, annot_fp))
        break
    out_dir = os.path.join(reg_dir, annot)
    os.makedirs(out_dir, exist_ok=True)
    out = os.path.join(out_dir, "ldsc_reg")
    if resume_prev_run and os.path.isfile(out+".results"):
        print("found exising file; skipped %s"%out)
        continue
    else:
        print(annot)
    subprocess.call([
        ldsc_bin,
        "--h2", sumstats,
        "--ref-ld-chr", "%s,%s"%(annot_fp,ldsc_baselineLD),
        "--frqfile-chr", ldsc_frqfile,
        "--w-ld-chr", ldsc_weight,
        "--overlap-annot",
        "--print-coefficients",
        "--print-delete-vals",
        "--out", out,
        "--chisq-max", "100000"]
        )


# 2. compute prop_h2g for each quantile of each continous annot
run_second_step = True
if run_second_step:
    print("2. compute prop_h2g for each quantile of each continous annot")
    for annot in annot_list:
        annot = annot.replace("|","--").replace("(","_").replace(")","_")
        annot_fp = os.path.join(par_dir, vep_dir, "label_wise_ldsc", "%sL2"%annot, "%sL2."%annot)
        out1 = os.path.join(reg_dir, annot, "quantile.q5.M")
        out2 = os.path.join(reg_dir, annot, "ldsc_reg.q5.txt")
        out0 = os.path.join(reg_dir, annot, "ldsc_reg")
        
        if resume_prev_run and os.path.isfile(out1):
            print("found existing file; skipped %s"%out1)
        else:
            print(out1)
            subprocess.call([
                    ldsc_quantile_M_bin,
                    "--ref-annot-chr", "%s,%s"%(annot_fp,ldsc_baselineLD),
                    "--frqfile-chr", ldsc_frqfile,
                    "--nb-quantile", "5",
                    "--maf", "0.05",
                    "--exclude0",
                    "--annot-header", annot,
                    "--out", out1]
                    )
        
        if resume_prev_run and os.path.isfile(out2):
            print("found existing file; skipped %s"%out2)
        else:
            print("compute h2g..")
            subprocess.call([
                ldsc_quantile_h2g_bin,
                out1,
                out0,
                out2
                ])

# 3. log done
print("3. log done")
