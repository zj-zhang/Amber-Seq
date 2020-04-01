"""read in the final stage training models, evaluation, and
plot their performances
ZZ
2020.2.10
"""

import os
import sys
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from src.nas_final.scorer import main as scorer_main
from src.utils import run_from_ipython
import argparse


def read_keras_new(fp):
    # updated keras should have the "metrics.log" 
    # that documents AUC-related metrics.
    dir_path = os.path.dirname(fp)
    base_metrics = read_performance_metrics(os.path.join(dir_path, "metrics.log"))
    with open(fp, "r") as f:
        for line in f:
            line = line.strip()
            if not line.startswith("Total time elapsed"):
                continue
            total_time = float(line.split(": ")[1])
    base_metrics['train_time'] = total_time
    return base_metrics


def read_performance_metrics(fp):
    record_order = ['val_auc', 'val_aupr', 'val_loss', 'test_auc', 'test_aupr', 'test_loss']
    metrics = {}
    with open(fp, "r") as f:
        i = 0
        for line in f:
            ele = line.strip().split()
            metrics[record_order[i]] = float(ele[2].strip('()[],'))
            i += 1
    return metrics


def read_keras_params(fp):
    with open(fp, "r") as f:
        for line in f:
            if not line.startswith("Trainable"):
                continue
            num = line.strip().split()[-1]
            num = int("".join(num.split(",")))
            return num


def read_batch_run(par_dir, add_baseline=False):
    #runs = [x for x in os.listdir(par_dir) if x.startswith("tmp_final")]
    runs = ['nas_final', 'nas_sample']
    vals_auc = []
    vals_aupr = []
    tests_auc = []
    tests_aupr = []
    category = []
    params = []
    train_time = []
    for r in runs:
       # use updated folder structures reader; zz 2020.2.19
        fp = "log.txt"
        metrics = read_keras_new(os.path.join(par_dir, r, fp))
        vals_auc.append(metrics['val_auc'])
        vals_aupr.append(metrics['val_aupr'])
        tests_auc.append(metrics['test_auc'])
        tests_aupr.append(metrics['test_aupr'])
        if r.endswith(".noSearch"):
            category.append("noSearch")
        elif r.endswith(".Random"):
            category.append("Random")
        else:
            category.append("Search")
        param_num = read_keras_params(os.path.join(par_dir, r, "keras_model_summary.txt"))
        params.append(param_num)
        train_time.append(metrics['train_time'])

    df = pd.DataFrame({
        "runs": runs,
        "vals_auc": vals_auc,
        "vals_aupr": vals_aupr,
        "tests_auc": tests_auc,
        "tests_aupr": tests_aupr,
        "category": category,
        "params": params,
        "train_time": train_time
        })
    if add_baseline:
        df = df.append(pd.Series({'runs':'baseline', 'vals_aupr':0.402, 'tests_aupr':0.338, 'vals_auc':0.937, 'tests_auc':0.931, 'params':52843120, 'category':'baseline', 'train_time':2073600.0}, name="baseline"))
    return df


def plot_batch_run_by_params(df, savefn=None):
    plt.close('all')
    fig, axs = plt.subplots(2, 2, figsize=(10,10))

    sns.lineplot(x="params", y="tests_aupr", hue="category", style="category", data=df, ax=axs[0,0])
    axs[0][0].set_title("Testing AUPR")

    sns.lineplot(x="params", y="vals_aupr", hue="category", style="category", data=df, ax=axs[0,1])
    axs[0][1].set_title("Validation AUPR")

    sns.lineplot(x="params", y="tests_auc", hue="category", style="category", data=df, ax=axs[1,0])
    axs[1][0].set_title("Testing AUC")

    sns.lineplot(x="params", y="vals_auc", hue="category", style="category", data=df, ax=axs[1,1])
    axs[1][1].set_title("Validation AUC")

    if savefn:
        plt.savefig(savefn)
    else:
        return axs


def plot_batch_run_by_order(df, savefn=None):
    df['rank'] = df.groupby("category")['tests_aupr'].rank()
    fig, ax = plt.subplots(1,1, figsize=(5,5))
    sns.lineplot(x="rank", y="tests_aupr", hue="category", style="category", data=df, ax=ax)
    if savefn:
        plt.savefig(savefn)
    else:
        return ax


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--par-dir", type=str, help="parent folder path")
    parser.add_argument("--od", type=str, help="output folder")
    args = parser.parse_args()

    od = args.od
    os.makedirs(od, exist_ok=True)
    df = read_batch_run(args.par_dir)
    plot_batch_run_by_order(df, os.path.join(od, "raw.png"))
    df.to_csv(os.path.join(od, "raw.tsv"), sep="\t", index=False)

    for i in range(df.shape[0]):
        this = df.iloc[i]['runs']
        scorer_main(model_fp=os.path.join(args.par_dir, this, "bestmodel.h5"), fn=os.path.join(args.od, "%s.scorer_output.txt"%this))



if __name__ == "__main__":
    if run_from_ipython():
        # running interactively
        print("running interactively")
        pass
    else:
        print("running in cmd")
        main()
