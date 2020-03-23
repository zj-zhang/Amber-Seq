# -*- coding: utf-8 -*-

import tensorflow as tf
import numpy as np
from keras.models import load_model
from .read_data import read_val_data, read_train_data, read_test_data, read_label_annot

from BioNAS.Preprocessing.pipeline import CnnFeatureModel
from BioNAS.Preprocessing.cluster import split_multitask_labels

from BioNAS.KFunctions.MotifKnowledgeFunc import MotifSaliency, make_output_annot

from BioNAS.Controller.reward import Loss_Reward, LossAucReward
from BioNAS.Controller.model_space import State, ModelSpace
from BioNAS.Controller.common_ops import batchify

from pkg_resources import resource_filename
from keras.callbacks import Callback
from keras.optimizers import Adam, SGD

DeepSEA_opt = SGD(lr=0.1, momentum=0.9, decay=1e-6, nesterov=False)


def read_controller_train_history(fn, last_only=None, is_randomized=False):
    import pandas as pd
    d = pd.read_csv(fn, sep=",")
    if last_only is not None:
        d = d.loc[ d.iloc[:,0] >= (d.shape[0] - last_only)]
    best_auc = 0
    best_arc = None
    if is_randomized:
        real_auc = [float(x.split(",")[1][:-1]) for x in d.iloc[:,1]]
        d.iloc[:,2] = real_auc
    for i in range(d.shape[0]):
        index, auc = d.iloc[i][0], d.iloc[i][2]
        arc = d.iloc[i][3:].to_list()
        if  auc > best_auc:
            best_arc = arc
            best_auc = auc
    return best_arc, best_auc


def prepare_base_model_and_data(session, feature_assign_fn=None, val_only=False, trainable=False):
    label_annot = read_label_annot()
    feature_assign_fn = feature_assign_fn or "../prep/channels_kmean_assign.npy" # default
    cat_list = ['TF', 'Pol', 'DNase', 'Histone']
    split_vec = [np.where(label_annot.category==x)[0] for x in cat_list]
    output_nodes = [State('Dense', activation="sigmoid", units=len(split_vec[i]), name="output_%s"%cat_list[i])
                    for i in range(len(split_vec))]

    x_val, y_val_ = read_val_data()
    y_val = split_multitask_labels(y_val_, split_vec)
    y_shapes = [(None,)+t.shape[1:] for t in y_val]

    if val_only:
        data_pack = [(x_val, y_val), (x_val, y_val)]
    else:
        x_train, y_train_ = read_train_data()
        y_train = split_multitask_labels(y_train_, split_vec)

        data_pack = [(x_train, y_train), (x_val, y_val)]

    feature_model_name = 'CnnFeatureModel'
    with tf.variable_scope(feature_model_name):
        cnn_model = load_model('deepsea_keras/data/deepsea_cnn.h5')
        feature_model = CnnFeatureModel(cnn_model, session=session, feature_assign_fn=feature_assign_fn,
                                        trainable=trainable,
                                        name=feature_model_name)
    return feature_model, output_nodes, data_pack


def get_model_space(num_layers):
    state_space = ModelSpace()
    for i in range(num_layers):
        state_space.add_layer(i, [
            State('Dense', units=50, activation='relu'),
            State('Dense', units=100, activation='relu'),
            State('Dense', units=250, activation='relu'),
            State('Dense', units=500, activation='relu')
        ])
    return state_space


def get_reward_fn(session, lbd=1, loss_c=None, knowledge_c=None, verbose=0):

    label_annot = read_label_annot()
    cat_list = ['TF', 'Pol', 'DNase', 'Histone']
    output_annot = make_output_annot(label_annot, cat_list)
    msk = MotifSaliency(output_annot, session,
                        pos_prop=0.9,
                        neg_prop=0.1,
                        batch_size=75,
                        index_to_letter={0: 'A', 1: 'G', 2: 'C', 3: 'T'},
                        filter_motif=True,
                        normalize_to_size=20,
                        verbose=verbose)
    #motif_file = '../../BioNAS/resources/rbp_motif/encode_motifs.txt.gz'
    motif_file = resource_filename('BioNAS.resources', 'rbp_motif/encode_motifs.txt.gz')
    msk.knowledge_encoder(motif_file=motif_file)

    reward_fn = LossAucReward(method='aupr', knowledge_function=msk, Lambda=lbd, loss_c=loss_c, knowledge_c=knowledge_c)

    return reward_fn


def get_test_data():
    x_test, y_test_ = read_test_data()
    label_annot = read_label_annot()
    cat_list = ['TF', 'Pol', 'DNase', 'Histone']

    split_vec = [np.where(label_annot.category==x)[0] for x in cat_list]
    y_test = split_multitask_labels(y_test_, split_vec)
    test_data = (x_test, y_test)
    return test_data


class ROCCallback(Callback):
    def __init__(self, validation_data, method='aupr', verbose=0, use_secondary=False):
        super().__init__()
        self.validation_data_ = validation_data
        self.scorer = LossAucReward(method=method)
        if method == 'aupr':
            self.secondary_scorer = LossAucReward(method='auc')
        else:
            self.secondary_scorer = LossAucReward(method='aupr')
        self.verbose = verbose
        self.use_secondary = use_secondary

    def on_epoch_end(self, epoch, logs=None):
        if logs is not None:
            auc = self.scorer(self.model, self.validation_data_)[0]
            logs['val_auc'] = auc
            if self.use_secondary:
                auc2 = self.secondary_scorer(self.model, self.validation_data_)[0]
                logs['val_auc2'] = auc2
            if self.verbose:
                if self.use_secondary:
                    print("Epoch %i, auc=%.5f, auc2=%.5f"%(epoch, auc, auc2))
                else:
                    print("Epoch %i, auc=%.5f"%(epoch, auc))


import datetime
class TimeBudgetCallback(Callback):
    def __init__(self, time_budget, verbose=0):
        super().__init__()
        self.time_budget = time_budget
        self.budget_in_secs = None
        self.train_start_time = None
        self.train_end_time = None
        self.consumed_time = None
        self.verbose = verbose
        self.budget_in_secs = sum(x * int(t) for x, t in zip([3600, 60, 1], self.time_budget.split(":")))

    def on_train_begin(self, logs=None):
        self.train_start_time = datetime.datetime.now()
        if self.verbose:
            print("Time budget set to: %s"%self.time_budget)

    def on_epoch_end(self, epoch, logs=None):
        self.consumed_time = (datetime.datetime.now() - self.train_start_time).total_seconds()
        if self.verbose:
            print("used time: %.2f %%"%(self.consumed_time / self.budget_in_secs * 100) )
        if self.consumed_time >= self.budget_in_secs:
            print("training ceased because run out of time budget")
            self.model.stop_training = True
    
    def on_train_end(self, logs=None):
        print("Total time elapsed: %s"%self.consumed_time)


def make_generator(x, y, batch_size):
    while True:
        try:
            g = batchify(x, y, batch_size)
            for x_, y_ in g:
                yield x_, y_
        except StopIteration:
            pass
