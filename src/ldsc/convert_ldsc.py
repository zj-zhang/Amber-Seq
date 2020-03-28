from BioNAS.Interpret.heritability import read_vep, convert_to_ldsc_annot, split_labels_to_folders
import sys
import os
import shutil
import subprocess
import tempfile
import gc

import argparse

parser = argparse.ArgumentParser(description="Step 1 of LDSC: convert annotation to L2")
parser.add_argument("--vep-dir", type=str, help="filepath to vep folder; will automatically find files by suffix")
parser.add_argument("--out-dir", type=str, help="filepath to folder for outputing L2 files")
parser.add_argument("--chroms", type=str, help="chroms (w/o chr) for the current task")
parser.add_argument("--label", type=str, help="filepath to the multi-task label annotations")
parser.add_argument("--ldsc-bfile", type=str, help="filepath to LDSC bfile from plink")
parser.add_argument("--ldsc-snps", type=str, help="filepath to LDSC SNPs list (e.g. listHM3.txt)")
parser.add_argument("--ldsc-baselineLD", type=str, help="filepath to LDSC baselineLD folder")
parser.add_argument("--ldsc-bin", type=str, help="filepath to executable (chmod +x) LDSC python script")

args = parser.parse_args()

# parsed argument
vep_dir = args.vep_dir
chroms_part = args.chroms
vep_dir = os.path.realpath(vep_dir)
out_dir = args.out_dir

# script filepath constants
#label_fp = "/mnt/ceph/users/zzhang/DeepSEA/ldsc_resources/labels_annot/total_label_idx.csv"
label_fp = args.label
#ldsc_bfile = "./ldsc_resources/1000G_EUR_QC/1000G_EUR_Phase3_plink/1000G.EUR.QC" # for baselineLD v2.2 and baseline v1.2
ldsc_bfile = args.ldsc_bfile
#ldsc_snps = "./ldsc_resources/HapMap3_SNPs/listHM3.txt"
baselineLD_dir = args.ldsc_baselineLD
ldsc_snps = args.ldsc_snps
ldsc_temp_dir_obj = tempfile.TemporaryDirectory() # node-local temp fs
ldsc_temp_dir = ldsc_temp_dir_obj.name

# first step: match abs_diffs to baselineLD SNPs
print("-"*40+' 1 '+"-"*40)
print("converting to ldsc-annot..")
convert_to_ldsc_annot(
        vep_dir=vep_dir,
        label_fp=label_fp,
        #baselineLD_dir="/mnt/home/zzhang/ceph/human_SNP/baselineLD/",
        baselineLD_dir=baselineLD_dir,
        output_dir=vep_dir,
        chroms_part=chroms_part,
        use_temp=ldsc_temp_dir
        )
gc.collect()

# second step: call ldsc to compute ld-score
#ldsc_bin = "./ldsc_resources/ldsc/ldsc.py"
ldsc_bin = args.ldsc_bin
print("-"*40+' 2 '+"-"*40)
print("run ldsc--l2..")
for chrom in chroms_part.split(','):
    print(chrom)
    subprocess.call(
            [
                ldsc_bin,
                "--l2",
                "--bfile",
                "%s.%s"%(ldsc_bfile, chrom),
                "--ld-wind-cm" ,
                "1.0",
                "--print-snps",
                "%s"%ldsc_snps,
                "--annot",
                "%s"%( os.path.join(ldsc_temp_dir, "joint.{chrom}.annot".format(chrom=chrom) )),
                "--out",
                "%s"%( os.path.join(ldsc_temp_dir, "joint.{chrom}".format(chrom=chrom) ))
            ]
            )

# third step: split into label-wise folders
print("-"*40+' 3 '+"-"*40)
print("split into output dir..")
os.makedirs(os.path.join(out_dir, "label_wise_ldsc"), exist_ok=True)
for chrom in chroms_part.split(','):
    print(chrom)
    split_labels_to_folders(
            l2_prefix=os.path.join(ldsc_temp_dir, "joint.{chrom}".format(chrom=chrom) ),
            output_dir=os.path.join(out_dir, "label_wise_ldsc")
            )
    shutil.move(os.path.join(ldsc_temp_dir, "joint.%s.log"%chrom), os.path.join(vep_dir, "joint.%s.log"%chrom))

# forth step: clean up
print("-"*40+' 4 '+"-"*40)
remove_prefix = tuple(['joint.%s'%x for x in chroms_part.split(',')])
remove_files = [os.path.join(ldsc_temp_dir, x) for x in os.listdir(ldsc_temp_dir) if x.startswith(remove_prefix)]
print("remove following temp files:")
for f in remove_files:
    try:
        os.remove(f)
        print(f)
    except Exception as e:
        print("Failed for %s because of %s"%(f, e))
