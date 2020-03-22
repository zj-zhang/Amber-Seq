import tensorflow as tf
import numpy as np
import keras.backend as K
import os
from BioNAS.Controller.model_builder import EnasCnnModelBuilder
from BioNAS.Controller.general_controller import GeneralController
from BioNAS.Controller.model_space import State, ModelSpace
from BioNAS.Controller.common_ops import count_model_params

from deepsea_keras.read_data import read_val_data, read_train_data, read_test_data, read_label_annot
from keras.callbacks import ModelCheckpoint, EarlyStopping, LearningRateScheduler
from common_func_msk import ROCCallback, TimeBudgetCallback
from keras.optimizers import Adam, SGD

from BioNAS.Controller.manager import EnasManager
from BioNAS.Controller.train_env import EnasTrainEnv, setup_logger
from BioNAS.Controller.reward import LossAucReward
import logging
from search_conv import get_controller

import argparse
import pickle
import sys
from convert_keras import read_controller_train_history
import pandas as pd


def lr_schedule(epoch):
    """Learning Rate Schedule

    Learning rate is scheduled to be reduced after 80, 120, 160, 180 epochs.
    Called automatically every epoch as part of callbacks during training.

    # Arguments
        epoch (int): The number of epochs

    # Returns
        lr (float32): learning rate
    """
    lr = 1e-3
    if epoch > 80:
        lr *= 0.5e-3
    elif epoch > 60:
        lr *= 1e-3
    elif epoch > 40:
        lr *= 1e-2
    elif epoch > 20:
        lr *= 1e-1
    print('Learning rate: ', lr)
    return lr


def random_sample_controller(controller_config_fp, model_space):
    config = pickle.load(open(controller_config_fp, "rb"))
    skip_target = config['skip_target']
    arc_seq = []
    num_layers = len(model_space)
    num_skips = sum([layer_id for layer_id in range(num_layers)])
    skips = np.zeros( num_skips, dtype=int)
    skips_idx = np.random.choice(np.arange(num_skips), size=int(num_skips*skip_target), replace=False)
    skips[skips_idx] = 1
    skips = skips.tolist()
    #print(skips)
    print(np.mean(skips))
    for layer_id in range(num_layers):
        count = np.random.choice(np.arange(len(model_space[layer_id])), size=1, replace=False).tolist()
        if layer_id > 0:
            skip = [skips.pop(0) for _ in range(layer_id)]
            count.extend(skip)
        arc_seq.extend(count)
    return arc_seq


def main(sd, od, wsf, gap=False, disable_controller=False, verbose=1):
    config = tf.ConfigProto()
    config.gpu_options.allow_growth = True
    session = tf.Session(config=config)
    K.set_session(session)
    model_space = pickle.load(open(os.path.join(sd, "model_space.pkl"), "rb"))
    for l in model_space:
        for x in l:
            x.Layer_attributes['filters'] *= wsf

    if disable_controller:
        print("random sampled from controller")
        arc_seq = random_sample_controller(os.path.join(sd, "controller_config.pkl"), model_space)
        print(arc_seq)
    else:
        is_randomized = True if sd.strip("/").endswith("noSearch") else False
        arc_seq, best_auc = read_controller_train_history(
                os.path.join(sd, "train_history.csv"),
                last_only=100,
                is_randomized=is_randomized
                )
        print(arc_seq)
        print("best_auc=",best_auc)


    wd = od

    input_node = State('input', shape=(1000, 4), name="input", dtype=tf.float32)
    output_node = State('dense', units=919, activation='sigmoid')
    model_compile_dict = {
        'loss': 'binary_crossentropy',
        'optimizer': 'adam',
        #'optimizer': SGD(lr=0.1, momentum=0.9, decay=1e-6, nesterov=True),
        #'metrics': []
    }

    child_batch_size = 1000
    model_fn = EnasCnnModelBuilder(
        dag_func='EnasConv1dDAG',
        batch_size=child_batch_size,
        session=session,
        model_space=model_space,
        inputs_op=[input_node],
        output_op=[output_node],
        num_layers=len(model_space),
        l1_reg=1e-8,
        l2_reg=5e-7,
        model_compile_dict=model_compile_dict,
        controller=None,
        dag_kwargs={
            'train_fixed_arc': True,
            'fixed_arc': arc_seq,
            'stem_config':{
                    'flatten_op': 'global_avg_pool' if gap else 'flatten',
                    'fc_units': 925
                }
            }
        )

    vars_list = [v for v in tf.trainable_variables() if v.name.startswith(model_fn.dag.name)]
    # remove optimizer related vars (e.g. momentum, rms)
    vars_list = [v for v in vars_list if not v.name.startswith("%s/compile"%model_fn.dag.name)]
    print("total model params: ", count_model_params(vars_list))
    with open(os.path.join(od,"tensor_vars.txt"), "w") as f:
        for v in vars_list:
            f.write("%s\t%i\n"%(v.name, int(np.prod(v.shape).value) ))

    val_data = read_val_data()
    #train_data = val_data
    #test_data = val_data
    train_data = read_train_data()
    test_data = read_test_data()

    model_weight_fp = os.path.join(wd, 'best_model.h5')

    roc_callback = ROCCallback(validation_data=val_data, method='auc', use_secondary=True, verbose=1)

    model_checkpointer = ModelCheckpoint(
        model_weight_fp,
        monitor='val_auc',
        mode='max',
        save_best_only=True,
        save_weights_only=True,
        verbose=1)

    early_stopper = EarlyStopping(
        monitor='val_auc',
        mode='max',
        patience=20,
        verbose=1)

    lr_scheduler = LearningRateScheduler(lr_schedule)

    time_budget = TimeBudgetCallback(time_budget="72:00:00", verbose=1)

    epochs = 100
    model = model_fn()
    try:
        model.fit(train_data[0], train_data[1],
                  epochs=epochs,
                  batch_size=child_batch_size,
                  validation_data=val_data,
                  callbacks=[roc_callback, model_checkpointer, early_stopper, lr_scheduler, time_budget],
                  verbose=verbose)
    except KeyboardInterrupt:
        print("user interruptted; computing test evaluations..")

    model.load_weights(model_weight_fp)
    print("training complete; performing evaluation..")

    def logging(logger, line):
        logger.write(line + '\n')
        logger.flush()
        print(line)

    reward_fn = LossAucReward(method='auc')
    reward_fn2 = LossAucReward(method='aupr')
    with open('%s/metrics.log'%wd, 'w') as logger:
        line = "val_data auc: %s" % str(reward_fn(model, val_data))
        logging(logger, line)

        line = "val_data aupr: %s" % str(reward_fn2(model, val_data))
        logging(logger, line)

        line = "val_data loss: %s" % str(model.evaluate(val_data[0], val_data[1]))
        logging(logger, line)

        line = "test_data: %s" % str(reward_fn(model, test_data))
        logging(logger, line)

        line = "test_data aupr: %s" % str(reward_fn2(model, test_data))
        logging(logger, line)

        line = "test_data loss: %s" % str(model.evaluate(test_data[0], test_data[1]))
        logging(logger, line)





if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Command-line arguments for BioNAS final train')
    
    parser.add_argument("--od", dest='od', type=str,
                        help='output dir')
    
    parser.add_argument('--sd', dest='sd', type=str,
                        help='search dir')
    
    parser.add_argument('--wsf', dest='width_scale_factor', default=1,
                        type=int,
                        help='width scale fatcor: multiplier for num of filters in each layer')

    parser.add_argument('--GAP', dest='gap', default=False,
                        action="store_true",
                        help='use GlobalAveragePool instead of Flatten')
    
    parser.add_argument('--disable-controller', dest='disable_controller',
                        default=False, action="store_true",
                        help='disable controller training')
    
    parser.add_argument('--verbose', dest='verbose', default=1,
                        type=int,
                        help='verbose mode')

    args = parser.parse_args()
    if args.od:
        main(sd=args.sd, od=args.od, wsf=args.width_scale_factor, gap=args.gap, disable_controller=args.disable_controller, verbose=args.verbose )
