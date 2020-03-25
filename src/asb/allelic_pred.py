import h5py
import sys
import os
import gc

import pandas as pd
import numpy as np
# use R for plotting.. 2020.2.24
#import matplotlib.pyplot as plt

import argparse

parser = argparse.ArgumentParser(description="Map Allelic-specific events from h5py")
parser.add_argument("--vep-prefix", type=str, help="Prefix for variant effect prediction (.h5) filepath")
parser.add_argument("--ref-pkl", type=str, help="reference mapping from ASB events to col/row index (.pkl)")
parser.add_argument("--out", type=str, help="filepath for output TSV file")

args = parser.parse_args()

asb_df = pd.read_pickle(args.ref_pkl)


pred_dict = {'Ref_pred': '%s.ref_predictions.h5'%args.vep_prefix,
        'Alt_pred':   '%s.alt_predictions.h5'%args.vep_prefix}

for k,v in pred_dict.items():
    print(k)
    pred_fp = v
    store = h5py.File(pred_fp, 'r')
    pred_all = store['data'].value
   
    pred = [pred_all[int(asb_df['row_idx'].iloc[i]), int(asb_df['col_idx'].iloc[i])] for i in range(asb_df.shape[0]) ]
    asb_df[k] = pred

    
    del pred_all
    store.close()
    gc.collect()


asb_df.to_csv(args.out, index=False, sep="\t")
