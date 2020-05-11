"""
Convert a particular arc_seq from a model_space to a Keras Model
and train
ZZJ
2020.1.10
"""
#from __future__ import absolute_import

import os
import argparse
import pickle
import numpy as np
from src.utils.model_spacer import get_model_space
from keras.models import Model
from keras.layers import Conv1D, Activation, BatchNormalization, Input, \
    MaxPooling1D, AveragePooling1D, Dropout, GlobalAveragePooling1D, Lambda, Add, Dense, Flatten
from keras import regularizers
from keras.utils import plot_model
from BioNAS.Controller.model_space import State, ModelSpace
from BioNAS.utils.plots import plot_training_history
from BioNAS.Controller.reward import LossAucReward
from contextlib import redirect_stdout
import json

# for training
from src.utils.read_data import read_val_data, read_train_data, read_test_data, read_label_annot
from keras.callbacks import ModelCheckpoint, EarlyStopping, LearningRateScheduler, ReduceLROnPlateau
from src.utils.common_func_msk import ROCCallback, TimeBudgetCallback, read_controller_train_history
from keras.optimizers import Adam, SGD
from keras.constraints import max_norm
from keras.utils import multi_gpu_model


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


def factorized_reduction_layer(input, out_filter, name, reduction_factor=4):
    x = Conv1D(out_filter,
               kernel_size=1,
               strides=1,
               kernel_initializer='he_normal',
               use_bias=False,
               padding="same",
               name=name
               )(input)
    x = MaxPooling1D(pool_size=reduction_factor, strides=reduction_factor, padding="same")(x)
    return x


def res_layer(layer, width_scale_factor, inputs, l2_reg=5e-7, name="layer", add_conv1_under_pool=True):
    if layer.Layer_type == 'conv1d':
        activation = layer.Layer_attributes['activation']
        num_filters = width_scale_factor * layer.Layer_attributes['filters']
        kernel_size = layer.Layer_attributes['kernel_size']
        if 'dilation' in layer.Layer_attributes:
            dilation = layer.Layer_attributes['dilation']
        else:
            dilation = 1
        x = Conv1D(num_filters,
                   kernel_size=kernel_size,
                   strides=1,
                   padding='same',
                   kernel_initializer='he_normal',
                   kernel_regularizer=regularizers.l2(l2_reg),
                   kernel_constraint=max_norm(0.9),
                   use_bias=False,
                   name="%s_conv"%name if dilation==1 else "%s_conv_d%i"%(name, dilation),
                   dilation_rate=dilation
                   )(inputs)
        x = BatchNormalization(name="%s_bn"%name)(x)
        if activation in ("None", "linear"):
            pass
        elif activation in ("relu", "sigmoid", "tanh", "softmax", "elu"):
            x = Activation(activation, name="%s_%s"%(name, activation))(x)
        elif activation == "leaky_relu":
            from keras.layers import LeakyReLU
            x = LeakyReLU(alpha=0.2, name="%s_%s"%(name, activation) )(x)
        else:
            raise Exception("Unknown activation: %s"%activation)
    elif layer.Layer_type == 'maxpool1d' or layer.Layer_type == 'avgpool1d':
        num_filters = width_scale_factor * layer.Layer_attributes['filters']
        pool_size = layer.Layer_attributes['pool_size']
        if add_conv1_under_pool:
            x = Conv1D(num_filters,
                       kernel_size=1,
                       strides=1,
                       padding='same',
                       kernel_initializer='he_normal',
                       use_bias=False,
                       name="%s_maxpool_conv"%name
                       )(inputs)
            x = BatchNormalization(name="%s_bn"%name)(x)
            x = Activation("relu", name="%s_relu"%name)(x)
        else:
            x = inputs

        if layer.Layer_type == 'maxpool1d':
            x = MaxPooling1D(pool_size=pool_size, strides=1, padding='same', name="%s_maxpool"%name)(x)
        elif layer.Layer_type == 'avgpool1d':
            x = AveragePooling1D(pool_size=pool_size, strides=1, padding='same', name="%s_avgpool"%name)(x)
        else:
            raise Exception("Unknown pool: %s"%layer.Layer_type)

    elif layer.Layer_type == 'identity':
        x = Lambda(lambda t: t, name="%s_id"%name)(inputs)
        #x = inputs
    else:
        raise Exception("Unknown type: %s"%layer.Layer_type)
    return x


def get_out_filters(model_space):
    out_filters = []
    pool_layers = []
    for layer_id in range(len(model_space)):
        layer = model_space[layer_id]
        this_out_filters = [l.Layer_attributes['filters'] for l in layer]
        assert len(set(this_out_filters)) == 1, "EnasConv1dDAG only supports one identical number of filters per layer," \
                                                "but found %i in layer %s" % (len(set(this_out_filters)), layer)
        if len(out_filters) and this_out_filters[0] != out_filters[-1]:
            pool_layers.append(layer_id - 1)

        out_filters.append(this_out_filters[0])
    print(out_filters)
    print(pool_layers)
    return out_filters, pool_layers


def convert(arc_seq, model_space, width_scale_factor=3, dropout_rate=0.1, model_plot_fn=None, pooling_mode='Flatten', add_conv1_under_pool=True):
    out_filters, pool_layers = get_out_filters(model_space)

    input = Input(shape=(1000, 4), name="input")
    # this is assuming all choices have the same out_filters
    stem_conv = State('conv1d', kernel_size=8, filters=out_filters[0], activation="linear")
    x = res_layer(stem_conv, width_scale_factor, input, name="stem_conv", add_conv1_under_pool=add_conv1_under_pool)

    start_idx = 0
    layers = []
    for layer_id in range(len(model_space)):
        print("start_idx=%i, layer id=%i, out_filters=%i x %i" % (start_idx, layer_id, out_filters[layer_id], width_scale_factor) )
        count = arc_seq[start_idx]
        this_layer = model_space[layer_id][count]
        print(this_layer)
        if layer_id == 0:
            x = res_layer(this_layer, width_scale_factor, x, name="L%i"%layer_id, add_conv1_under_pool=add_conv1_under_pool)
        else:
            x = res_layer(this_layer, width_scale_factor, layers[-1], name="L%i"%layer_id, add_conv1_under_pool=add_conv1_under_pool)

        if layer_id > 0:
            skip = arc_seq[start_idx+1 : start_idx+layer_id+1]
            skip_layers = [layers[i] for i in range(len(layers)) if skip[i]==1]
            print("skip=%s"%skip)
            if len(skip_layers):
                skip_layers.append(x)
                x = Add(name="L%i_resAdd"%layer_id)(skip_layers)
            x = BatchNormalization(name="L%i_resBn"%layer_id)(x)

        if dropout_rate != 0:
            x = Dropout(dropout_rate, name="L%i_dropout"%layer_id)(x)

        layers.append(x)
        if layer_id in pool_layers:
            pooled_layers = []
            for i, layer in enumerate(layers):
                pooled_layers.append(
                    factorized_reduction_layer(
                        layer,
                        out_filters[layer_id+1]*width_scale_factor,
                        name="pool_at_%i_from_%i"%(layer_id, i))
                )
            print("pooled@%i, %s"%(layer_id, pooled_layers))
            layers = pooled_layers

        start_idx += 1 + layer_id
        print('-'*80)

    # fully-connected layer
    if pooling_mode == 'GAP':
        x = GlobalAveragePooling1D()(x)
    elif pooling_mode == 'Flatten':
        x = Flatten()(x)
    else:
        raise Exception("Unknown pooling mode: %s"%pooling)
    x = Dropout(dropout_rate)(x)
    x = Dense(units=925, activation="relu")(x)

    output = Dense(units=919, activation="sigmoid", name="output")(x)

    model = Model(inputs=input, outputs=output)
    if model_plot_fn:
        plot_model(model, model_plot_fn)
    return model


def get_num_layers(arc_seq):
    start_idx = 0
    for layer_id in range(100):
        start_idx += 1 + layer_id
        if start_idx >= len(arc_seq):
            return layer_id+1
    raise Exception("Too many layers (>100)")


def main(sd, wd, final_config, gpus=1, pooling='Flatten', width_scale_factor=2, disable_controller=False, add_conv1_under_pool=True, verbose=1):
    # working directory
    if gpus > 1:
        use_multi_gpu = True
    else:
        use_multi_gpu = False
    if not os.path.isdir(wd):
        os.makedirs(wd)

    # load model space
    if os.path.isfile(os.path.join(sd, "model_space.pkl")):
        print("load pre-computed model_space")
        model_space = pickle.load(open(os.path.join(sd, "model_space.pkl"), "rb"))
    else:
        raise Exception("cannot find model_space.pkl file")

    # read the best architecture
    if disable_controller:
        print("random sampled from controller")
        arc_seq = random_sample_controller(os.path.join(sd, "controller_config.pkl"), model_space)
        print(arc_seq)
    else:
        is_randomized = True if sd.strip("/").endswith("noSearch") else False
        arc_seq, search_auc = read_controller_train_history(
                os.path.join(sd, "train_history.csv"),
                last_only=100,
                is_randomized=is_randomized
                )
        print("read search results from %s"% sd)
        print("best searched auc=%.5f"%search_auc)
        num_layers = get_num_layers(arc_seq)
        print("inferred %i layers" % num_layers)

    dropout_rate = final_config.pop("dropout_rate", 0.1)

    print("pooling mode is %s"%pooling)
    model = convert(arc_seq=arc_seq, 
                model_space=model_space, 
                width_scale_factor=width_scale_factor, 
                dropout_rate=dropout_rate,
                model_plot_fn=os.path.join(wd,"model.png"),
                pooling_mode=pooling,
                add_conv1_under_pool=add_conv1_under_pool
                )

    with open(os.path.join(wd, 'keras_model_summary.txt'), 'w') as f:
        with redirect_stdout(f):
            model.summary()

    # read datas
    val_data = read_val_data()
    #train_data = val_data
    # test_data = val_data
    train_data = read_train_data()

    # preps for training
    model_weight_fp = os.path.join(wd, "bestmodel.h5")
    roc_callback = ROCCallback(validation_data=val_data, method='auc', verbose=verbose, use_secondary=True)
    checkpointer = ModelCheckpoint(
        model_weight_fp,
        monitor='val_auc',
        mode='max',
        save_best_only=True,
        save_weights_only=False,
        verbose=verbose)
    earlystopper = EarlyStopping(
        monitor='val_auc',
        mode='max',
        patience=final_config.pop("early_stop_patience", 10),
        verbose=verbose)

    lr_scheduler = LearningRateScheduler(lr_schedule)

    lr_reducer = ReduceLROnPlateau(
        monitor='val_auc',
        mode='max',
        factor=np.sqrt(0.1),
        cooldown=0,
        patience=final_config.pop("reduce_lr_patience", 5),
        min_lr=0.5e-6,
        verbose=verbose
    )

    time_budget = TimeBudgetCallback(time_budget=final_config.pop("train_time_budget", "24:00:00"), verbose=1)

    # training
    if use_multi_gpu:
        model = multi_gpu_model(model, gpus=gpus)

    model.compile(
        loss='binary_crossentropy',
        optimizer='adam',
        metrics=['acc']
    )
    batch_size_multiplier = gpus
    try:
        hist = model.fit(
            train_data[0], train_data[1],
            epochs=final_config.pop("train_epochs", 200),
            batch_size=1000 * batch_size_multiplier,
            verbose=verbose,
            validation_data=val_data,
            callbacks=[roc_callback, checkpointer, earlystopper, lr_reducer, lr_scheduler, time_budget]
        )
        pickle.dump(hist, open(os.path.join(wd, "hist.pkl"), 'wb'), -1)
        plot_training_history(hist, wd)

    except KeyboardInterrupt:
        pass

    # testing
    test_data = read_test_data()
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


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Command-line arguments for Keras final train')
    
    parser.add_argument("--od", dest='od', type=str,
                        help='output dir')
    
    parser.add_argument('--sd', dest='sd', type=str,
                        help='search dir')
    
    parser.add_argument('--config', dest='config', type=str,
                        help='configuration filepath')
    
    parser.add_argument('--pooling', dest='pooling', type=str,
                        default='Flatten',
                        choices=['GAP', 'Flatten'],
                        help='Pooling mode for conversion from convolution to dense layers')
    
    parser.add_argument('--gpus', dest='gpus', default=1,
                        type=int,
                        help='number of gpus')

    parser.add_argument('--wsf', dest='width_scale_factor', default=1,
                        type=int,
                        help='width scale fatcor: multiplier for num of filters in each layer')

    parser.add_argument('--disable-controller', dest='disable_controller',
                        default=False, action="store_true",
                        help='disable controller training')

    parser.add_argument('--remove-conv1-under-pool', dest='remove_conv1_under_pool',
                        default=False, action="store_true",
                        help='Do Not add 1x1 conv before max/avg pooling')

    parser.add_argument('--verbose', dest='verbose', default=1,
                        type=int,
                        help='verbose mode')

    args = parser.parse_args()
    if args.od is not None:
        assert os.path.isfile(args.config), "Config file not found: %s"%args.config
        assert args.config.endswith("json"), "Must have a json file for configuration"
        final_config = json.load(open(args.config, "r"))
        main(
                wd=args.od,
                sd=args.sd,
                final_config=final_config,
                gpus=args.gpus,
                pooling=args.pooling,
                width_scale_factor=args.width_scale_factor,
                disable_controller=args.disable_controller,
                add_conv1_under_pool=(not args.remove_conv1_under_pool),
                verbose=args.verbose
                )
