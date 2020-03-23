import h5py
import sys
import os
import gc

import pandas as pd
import numpy as np
# use R for plotting.. 2020.2.24
#import matplotlib.pyplot as plt


model_dir = sys.argv[1]

# DEPRECATED: now use new pickle to read
#asb_df = pd.read_csv("./asb/asb_with_index.tsv", sep="\t")
#asb_df = asb_df.loc[(~np.isnan(asb_df['vep_row']) & (~np.isnan(asb_df['vep_col'])) )]

dgf_df = pd.read_pickle("./asb/DGF_pandas_table.pkl")
tf_df = pd.read_pickle("./asb/TF_pandas_table.pkl")

dfs = {"DGF": dgf_df, "TF":tf_df}

pred_dict = {'Ref_pred': 'hg19.baselineLD_All.20180423.ref_predictions.h5',
        'Alt_pred':   'hg19.baselineLD_All.20180423.alt_predictions.h5'}

for k,v in pred_dict.items():
    print(k)
    pred_fp = os.path.join(model_dir, 'vep_output', v)
    store = h5py.File(pred_fp, 'r')
    pred_all = store['data'].value
   
    for df_name, asb_df in dfs.items():
        print(df_name)
        pred = [pred_all[int(asb_df['row_idx'].iloc[i]), int(asb_df['col_idx'].iloc[i])] for i in range(asb_df.shape[0]) ]
        asb_df[k] = pred

    for df_name, asb_df in dfs.items():
        #asb_df = asb_df.drop(columns=['File name'])
        asb_df.to_csv(os.path.join(model_dir, "%s.Allelic_Specific.tsv"%df_name), index=False, sep="\t")
    
    del pred_all
    store.close()
    gc.collect()


