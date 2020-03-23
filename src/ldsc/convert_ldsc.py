from BioNAS.Interpret.heritability import read_vep, convert_to_ldsc_annot, split_labels_to_folders
import sys
import os
import shutil
import subprocess
import tempfile
import gc

# parsed argument
vep_dir = sys.argv[1]
chroms_part = sys.argv[2]
vep_dir = os.path.realpath(vep_dir)

# script filepath constants
label_fp = "/mnt/ceph/users/zzhang/DeepSEA/ldsc_resources/labels_annot/total_label_idx.csv"
ldsc_bfile = "./ldsc_resources/1000G_EUR_QC/1000G_EUR_Phase3_plink/1000G.EUR.QC" # for baselineLD v2.2 and baseline v1.2
#ldsc_bfile = "./ldsc_resources/1000G_EUR_QC/1000G_Phase1_plinkfiles/1000G_plinkfiles/1000G.mac5eur"  # for baseline v1.2
ldsc_snps = "./ldsc_resources/HapMap3_SNPs/listHM3.txt"
#ldsc_temp_dir_obj = tempfile.TemporaryDirectory() #"/scratch/"
#ldsc_temp_dir = ldsc_temp_dir_obj.name
ldsc_temp_dir = "./tmp"

# first step: match abs_diffs to baselineLD SNPs
print("-"*40+' 1 '+"-"*40)
print("converting to ldsc-annot..")
convert_to_ldsc_annot(
        vep_dir=vep_dir,
        label_fp=label_fp,
        #baselineLD_dir="/mnt/home/zzhang/ceph/human_SNP/baselineLD/",
        baselineLD_dir="/mnt/ceph/users/zzhang/human_SNP/1000G_Phase3_baseline_v1.1/baseline_v1.1/", # use baseline v1.2
        output_dir=vep_dir,
        chroms_part=chroms_part,
        use_temp=ldsc_temp_dir
        )
gc.collect()

# second step: call ldsc to compute ld-score
ldsc_bin = "./ldsc_resources/ldsc/ldsc.py"
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
os.makedirs(os.path.join(vep_dir, "label_wise_ldsc"), exist_ok=True)
for chrom in chroms_part.split(','):
    print(chrom)
    split_labels_to_folders(
            l2_prefix=os.path.join(ldsc_temp_dir, "joint.{chrom}".format(chrom=chrom) ),
            output_dir=os.path.join(vep_dir, "label_wise_ldsc")
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
