"""plotting functions for controller embedding etc
ZZ
2020.02.11
"""

import os
import sys
import tensorflow as tf
import pickle
import json
from src.nas_search.search_conv import get_controller
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from BioNAS.Controller.model_space import get_layer_shortname


def load_controller(controller, sd):
    assert os.path.isfile(os.path.join(sd, "model_space.pkl"))
    controller.load_weights(os.path.join(sd, "controller_weights.h5"))
    weights = controller.get_weights()
    weights_dict = {controller.weights[i].name: weights[i] for i in range(len(weights))}
    return weights_dict


def load_and_dump_nas_action_probs(sd, od):
    data = json.load(open(os.path.join(sd, "weight_data.json"),"r"))
    
    plt.close('all')
    fig, axes = plt.subplots(4, 4, figsize=(18, 18), dpi=200)
    #fig.tight_layout()
    fig.subplots_adjust(wspace=0.5, hspace=0.5)
    ax_set = [ axes[i,j] for i in range(4) for j in range(1,4)]
    ax_set = plot_action_weights_ax(data, ax_set)
    for i in range(4):
        handles, labels = axes[i,1].get_legend_handles_labels()
        axes[i,0].axis('off')
        axes[i,0].legend(handles, labels, loc='center')
    plt.savefig(os.path.join(od, "selection_prob_over_time.pdf"))
    
    sum_df = pd.DataFrame(columns=["layer_id", "layer_type", "prob"])
    for i in range(len(data)):
        d = {}
        for k in data['L%i'%i]['operation']:
            if k.startswith("conv"):
                new_k = k.split("_")
                new_k.pop(1)
                new_k = "_".join(new_k)
            else:
                new_k = k
            d[new_k] = data["L%i"%i]["operation"][k][0]
            sum_df = sum_df.append( {"layer_id": i, "layer_type": new_k, "prob": d[new_k][-1]}, ignore_index=True  )
        
        df = pd.DataFrame.from_dict(d)
        df.to_csv(os.path.join(od, "L%i.prob_over_time.tsv"%i), sep="\t", index=False)
    
    sum_df.to_csv(os.path.join(od, "selection_prob_at_convergence.tsv"), sep="\t", index=False)
    
    return data


def plot_action_weights_ax(df, ax_set):
    i = 0
    for layer in df:
        ax = ax_set[i]
        i += 1
        for type, data in df[layer]['operation'].items():
            d = np.array(data)
            avg = np.apply_along_axis(np.mean, 0, d)
            ax.plot(avg, label=type)
        ax.set_title(layer)
        ax.set_xlabel('Controller steps')
        ax.set_ylabel('Selection Prob')
    return ax_set


def princomp_embed(emb, model_space, pc_indices=None, savefn=None, save_data=False):
    assert emb.shape[0] == len(model_space[0])
    def label_point(x, y, val, ax):
        a = pd.DataFrame({'x': x, 'y': y, 'val': val})
        for i, point in a.iterrows():
            ax.text(point['x']+.02, point['y'], str(point['val']))

    if pc_indices is None:
        pc_indices = [0, 1]
    pca = PCA(n_components=4)
    pc = pca.fit_transform(emb)
    plt.close()
    fig, ax = plt.subplots(1,1)
    sns.scatterplot(x=pc[:,pc_indices[0]], y=pc[:,pc_indices[1]], ax=ax)
    shortnames = []
    for x in model_space[0]:
        n = get_layer_shortname(x) 
        if n.startswith("conv"):
            n = n.split("_")
            n.pop(1)
            n = "_".join(n)
        shortnames.append(n)
    layer_types = [x.Layer_type for x in model_space[0]] 
    label_point(pc[:,pc_indices[0]], pc[:,pc_indices[1]], val=shortnames, ax=ax)
    xlab = 'PC%i, %.2f' % (pc_indices[0]+1, pca.explained_variance_ratio_[pc_indices[0]]) 
    ax.set_xlabel(xlab)
    ylab =  'PC%i, %.2f' % (pc_indices[1]+1, pca.explained_variance_ratio_[pc_indices[1]]) 
    ax.set_ylabel(ylab)
    if savefn:
        plt.savefig(savefn)
        if save_data:
            datafn = os.path.join( os.path.dirname(savefn), "pca_%i-%i_data.txt"%(pc_indices[0]+1, pc_indices[1]+1) )
            df = pd.DataFrame({"labels": shortnames, xlab: pc[:, pc_indices[0]], ylab: pc[:,pc_indices[1]], 
                "layer_types": layer_types})
            df.to_csv(datafn, sep="\t", index=False)
    else:
        return ax


def tsne_embed(emb, model_space, savefn=None, random_state=None):
    assert emb.shape[0] == len(model_space[0])
    def label_point(x, y, val, ax):
        a = pd.DataFrame({'x': x, 'y': y, 'val': val})
        for i, point in a.iterrows():
            ax.text(point['x']+.02, point['y'], str(point['val']))

    dc = TSNE(n_components=2, n_iter=250, verbose=0, random_state=random_state).fit_transform(emb)
    dc_indices = [0, 1]
    plt.close()
    fig, ax = plt.subplots(1,1)
    sns.scatterplot(x=dc[:,dc_indices[0]], y=dc[:,dc_indices[1]], ax=ax)
    shortnames = [get_layer_shortname(x) for x in model_space[0]]
    label_point(dc[:,dc_indices[0]], dc[:,dc_indices[1]], val=shortnames, ax=ax)
    ax.set_xlabel('tSNE-1' )
    ax.set_ylabel('tSNE-2' )
    if savefn:
        plt.savefig(savefn)
    else:
        return ax


def main_batch_run(batch_run_dir):
    sds = [x for x in os.listdir(batch_run_dir) if x.startswith("tmp_search")]
    session = tf.Session()
    # WARNING:
    # assuming all individual runs under batch run will have the same model_space
    model_space = pickle.load(open(os.path.join(batch_run_dir, sds[0], "model_space.pkl"),"rb"))
    controller = get_controller(model_space, session)
    for sd in sds:
        print(sd)
        sd = os.path.join(batch_run_dir, sd)
        if not os.path.isfile(os.path.join(sd, "controller_weights.h5")):
            continue
        weights_dict = load_controller(controller, sd)
        emb = weights_dict["controller/create_weights/emb/layer_0/w_start:0"]
        princomp_embed(emb, model_space, pc_indices=[0,1], savefn=os.path.join(sd, "PCA_1-2_embed.png"))
        princomp_embed(emb, model_space, pc_indices=[1,2], savefn=os.path.join(sd, "PCA_2-3_embed.png"))
        # tsne is so random for only a few points.. it should be better for large clusters of data
        #tsne_embed(emb, model_space, savefn=os.path.join(sd, "tSNE_embed.png"), random_state=777)


def main(sd, od, controller_config):
    session = tf.Session()
    model_space = pickle.load(open(os.path.join(sd, "model_space.pkl"),"rb"))
    controller = get_controller(controller_config, model_space, session)
    weights_dict = load_controller(controller, sd)
    emb = weights_dict["controller/create_weights/emb/layer_0/w_start:0"]
    princomp_embed(emb, model_space, pc_indices=[0,1], savefn=os.path.join(od, "PCA_1-2_embed.png"), save_data=True)
    princomp_embed(emb, model_space, pc_indices=[1,2], savefn=os.path.join(od, "PCA_2-3_embed.png"), save_data=True)
    
    load_and_dump_nas_action_probs(sd, od)


if __name__ == "__main__":
    sd, od = sys.argv[1], sys.argv[2]
    if od != "None":
        print("analyze for a single run")
        controller_config = json.load(open(sys.argv[3], "r"))
        main(sd, od, controller_config)
    else:
        print("analyze for a batch run")
        batch_run_dir = sys.argv[1]
        if batch_run_dir:
             main_batch_run(batch_run_dir)
