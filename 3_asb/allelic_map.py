"""A belated script for mapping original Allele-specific binding events to the 
row- and column- indices in the hdf5 ref-alt prediction files
ZZ
2020.2.24
"""

import os
import sys
import h5py
import pandas as pd
import numpy as np
import pickle
from utils import deepsea_label_name_normalizer


def read_index_mapping():
    annot_df = pd.read_csv("./labels_annot/total_label_idx.csv", sep="\t")
    #vep_dict = pickle.load(open("./labels_annot/vep_dict.pkl", "rb"))
    vep_dict = pickle.load(open("./labels_annot/vep_dict_ref_alt.pkl", "rb"))
    return annot_df, vep_dict


def read_DGF_data(annot_df, vep_dict):
    data = pd.read_csv("asb/DeepSEA_ST3_allele_specific_DNase.csv", header=1)
    
    # alignment filename to cell-line info downloaded from ucsc-encode:
    # http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwDgf/
    aln_to_cell = {}
    with open("asb/encode_Dnase_name_mapping.txt", "r") as f:
        for line in f:
            ele = line.strip().split()
            if not ele[0].endswith("bam"):
                continue
            if ele[1].split("=")[1] == "revoked":
                continue
            aln = ele[0]
            cell = ele[7].split("=")[1].strip(";")
            aln_to_cell[aln] = cell
    cell_line_map = [aln_to_cell[ data['File name'][i] ] for i in range(data.shape[0]) ]
    data['new_cell_type'] = cell_line_map
    
    # now map to DeepSEA index
    annot_dict = { 
            x: annot_df.loc[(annot_df.cell==x) & (annot_df.target=='DNase')].label_idx.values
            for x in set(data['new_cell_type']) }
    col_idx = []
    row_idx = []
    for i in range(data.shape[0]):
        this_col = annot_dict[data['new_cell_type'][i]]
        assert len(this_col) == 1
        this_col = this_col[0]
        #row_id = (data['CHR'][i],str(data['POS'][i]))
        row_id = (data['CHR'][i], data['POS'][i], data['REF'][i], data['ALT'][i])
        if row_id in vep_dict:
            this_row = vep_dict[ row_id ]
        else:
            this_row = np.nan
        # there are multi-allelic variants that became multiple rows..
        #assert len(this_row) == 1, Exception("non-unique row index: %s, %s = %s"%(data['CHR'][i], str(data['POS'][i]), this_row))
        #this_row = this_row[0]
        col_idx.append(this_col)
        row_idx.append(this_row)

    data['col_idx'] = col_idx
    data['row_idx'] = row_idx
    #pickle.dump(data, open("./asb/DGF_pandas_table.pkl", "wb"))
    
    # drop some non-existing SNPs in hdf5.. might be in repeats
    new_data = data.dropna(subset=['row_idx'])
    print("removed %i non-existing SNPs from %i.."%(data.shape[0] - new_data.shape[0], data.shape[0]))
    data = new_data
    # drop multi-allelic SNPs
    #data = data.drop(axis=0, index=[data.index[i] for i in range(data.shape[0]) if len(data.iloc[i]['row_idx'])>1] )
    data['row_idx'] = data['row_idx'].astype('int')
    data.to_pickle("./asb/DGF_pandas_table.pkl")
    data.to_csv("./asb/DGF_pandas_table.tsv", sep="\t", index=False)
    return data


def read_TF_data(annot_df, vep_dict):
    data = pd.read_table("./asb/Frey-2018-Supp1.gz")
    data = data.dropna(subset=['deepsea_col'])
    data['deepsea_col'] = [deepsea_label_name_normalizer(x) for x in data['deepsea_col']]
    annot_dict = { '--'.join(x.split("--")[1:] ): int(x.split("--")[0]) for x in annot_df['label_name'] }
    old_annot_df =pd.read_table("./deepsea_keras/resources/label_names.original.txt", header=None)
    old_annot_dict = { deepsea_label_name_normalizer(old_annot_df.iloc[i,0]):i for i in range(old_annot_df.shape[0]) }

    col_idx = []
    row_idx = []
    for i in range(data.shape[0]):
        col_id = data.iloc[i]['deepsea_col']
        this_col = old_annot_dict[col_id]
        #row_id = (data.iloc[i]['chr'], str(data.iloc[i]['pos']) )
        row_id = (data.iloc[i]['chr'], data.iloc[i]['pos'], data.iloc[i]['ref'], data.iloc[i]['alt'] )
        if row_id in vep_dict:
            this_row = vep_dict[row_id]
        else:
            this_row = np.nan
        
        col_idx.append(this_col)
        row_idx.append(this_row)

    data['col_idx'] = col_idx
    data['row_idx'] = row_idx
    
    new_data = data.dropna(subset=['row_idx'])
    print("removed %i non-existing SNPs from %i.."%(data.shape[0] - new_data.shape[0], data.shape[0]))
    data = new_data
    data['row_idx'] = data['row_idx'].astype('int')
    data.to_pickle("./asb/TF_pandas_table.pkl")
    data.to_csv("./asb/TF_pandas_table.tsv", sep="\t", index=False)
    return data


def main():
    annot_df, vep_dict = read_index_mapping()
    d1 = read_DGF_data(annot_df, vep_dict)
    d2 = read_TF_data(annot_df, vep_dict)
    return d1, d2
