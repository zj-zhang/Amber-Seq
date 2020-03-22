import tensorflow as tf
import numpy as np
import keras.backend as K
import os
from BioNAS.Controller.model_builder import EnasCnnModelBuilder
from BioNAS.Controller.general_controller import GeneralController
from BioNAS.Controller.model_space import State, ModelSpace
from BioNAS.Controller.common_ops import count_model_params

from deepsea_keras.read_data import read_val_data, read_train_data, read_test_data, read_label_annot
from keras.callbacks import ModelCheckpoint, EarlyStopping
#from common_func_msk import ROCCallback
from keras.optimizers import Adam, SGD

from BioNAS.Controller.manager import EnasManager
from BioNAS.Controller.train_env import EnasTrainEnv, setup_logger
from BioNAS.Controller.reward import LossAucReward
from BioNAS.utils.plots import plot_controller_hidden_states
import logging
import pickle


#def get_reward_fn(session, lbd=1, loss_c=None, knowledge_c=None):
#    from BioNAS.KFunctions.MotifKnowledgeFunc import MotifSaliency, make_output_annot
#    from pkg_resources import resource_filename
#    label_annot = read_label_annot()
#    cat_list = [('TF', 'Pol', 'DNase', 'Histone')]
#    output_annot = make_output_annot(label_annot, cat_list)
#    print(output_annot[-1])
#    msk = MotifSaliency(output_annot, session,
#                        pos_prop=0.9,
#                        neg_prop=0.1,
#                        batch_size=100,
#                        index_to_letter={0: 'A', 1: 'G', 2: 'C', 3: 'T'},
#                        filter_motif=True,
#                        verbose=1)
#    motif_file = resource_filename('BioNAS.resources', 'rbp_motif/encode_motifs.txt.gz')
#    msk.knowledge_encoder(motif_file=motif_file)
#
#    reward_fn = LossAucReward(method='aupr', knowledge_function=msk, Lambda=lbd, loss_c=loss_c, knowledge_c=knowledge_c)
#
#    return reward_fn


def get_controller(model_space, session):
    with tf.device("/cpu:0"):
        #lr = 0.0 if disable_controller else 0.001
        lr = 0.001
        controller = GeneralController(
            model_space,
            session=session,
            share_embedding={i:0 for i in range(1, len(model_space))},
            with_skip_connection=True,
            with_input_blocks=False,
            num_input_blocks=1,
            skip_connection_unique_connection=False,
            skip_weight=1.0,
            skip_target=0.4,
            lstm_size=64,
            lstm_num_layers=1,
            kl_threshold=0.01,
            train_pi_iter=10,
            #optim_algo=SGD(lr=lr, momentum=True),
            optim_algo='adam',
            temperature=2.,
            lr_init=lr,
            lr_dec_rate=1.0,
            tanh_constant=1.5,
            buffer_size=1,  ## num of episodes saved
            batch_size=20
        )
        controller.buffer.rescale_advantage_by_reward = False
    return controller


def get_model_space(out_filters=64, num_layers=9, num_pool=4):
    state_space = ModelSpace()
    if num_pool == 4:
        expand_layers = [num_layers//4-1, num_layers//4*2-1, num_layers//4*3-1]
    elif num_pool == 3:
        expand_layers = [num_layers//3-1, num_layers//3*2-1]
    else:
        raise Exception("Unsupported pooling num: %i"%num_pool)
    for i in range(num_layers):
        state_space.add_layer(i, [
            State('conv1d', filters=out_filters, kernel_size=8, activation='relu'),
            State('conv1d', filters=out_filters, kernel_size=4, activation='relu'),
            State('conv1d', filters=out_filters, kernel_size=8, activation='relu', dilation=10),
            State('conv1d', filters=out_filters, kernel_size=4, activation='relu', dilation=10),
            # max/avg pool has underlying 1x1 conv
            State('maxpool1d', filters=out_filters, pool_size=4, strides=1),
            State('avgpool1d', filters=out_filters, pool_size=4, strides=1),
            State('identity', filters=out_filters),
      ])
        if i in expand_layers:
            out_filters *= 2

    return state_space


def get_controller_config(controller):
    d = controller.__dict__
    d = {k:d[k] for k in d  if type(d[k]) not in (tf.Tensor, tf.Variable, list, dict)}
    non_picklable_keys = ('optimizer', 'train_op', 'session', 'g_emb', 'w_attn_1', 'w_attn_2', 'v_attn', 'train_step', 'optim_algo' )
    d = {k:d[k] for k in d if not k in non_picklable_keys}
    print(d)
    return d


def main(wd, num_layers, gap, disable_controller, verbose):
    config = tf.ConfigProto()
    config.gpu_options.allow_growth = True
    session = tf.Session(config=config)
    K.set_session(session)
    out_filters = 32
    model_space = get_model_space(out_filters=out_filters, num_layers=num_layers, num_pool=4)
    controller = get_controller(model_space, session=session)
    pickle.dump(model_space, open(os.path.join(wd, "model_space.pkl"), "wb"))
    controller_config = get_controller_config(controller)
    pickle.dump(controller_config, open(os.path.join(wd, "controller_config.pkl"), "wb") )

    #reward_fn = LossAucReward(method='aupr')
    reward_fn = LossAucReward(method='auc')

    input_node = State('input', shape=(1000, 4), name="input", dtype=tf.float32)
    output_node = State('dense', units=919, activation='sigmoid')
    model_compile_dict = {
        'loss': 'binary_crossentropy',
        'optimizer': 'adam',
        #'optimizer': SGD(lr=0.1, momentum=0.9, decay=1e-6, nesterov=False),
        #'metrics': ['acc']
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
        controller=controller,
        dag_kwargs={
            'stem_config':{
                'flatten_op': 'global_avg_pool' if gap else 'flatten',
                'fc_units': 925
                },
            #'add_conv1_under_pool': False
            }
    )

    val_data = read_val_data()
    #train_data = val_data
    train_data = read_train_data()


    manager = EnasManager(
        train_data=train_data,
        validation_data=val_data,
        epochs=1,
        child_batchsize=child_batch_size,
        reward_fn=reward_fn,
        model_fn=model_fn,
        post_processing_fn='minimal',
        model_compile_dict=model_compile_dict,
        working_dir=wd,
        disable_controller=disable_controller,
        verbose=verbose
        )

    logger = setup_logger(wd, verbose_level=logging.DEBUG)
    vars_list = [v for v in tf.trainable_variables() if v.name.startswith(model_fn.dag.name)]
    # remove optimizer related vars (e.g. momentum, rms)
    vars_list = [v for v in vars_list if not v.name.startswith("%s/compile"%model_fn.dag.name)]
    logger.info("total model params: %i" % count_model_params(vars_list))

    with open(os.path.join(wd,"tensor_vars.txt"), "w") as f:
        for v in vars_list:
            f.write("%s\t%i\n"%(v.name, int(np.prod(v.shape).value) ))

    env = EnasTrainEnv(
        controller,
        manager,
        logger=logger,
        max_episode=300,
        max_step_per_ep=100,
        working_dir=wd,
        time_budget="72:00:00",
        with_input_blocks=False,
        with_skip_connection=True,
        child_train_steps=500,
        child_warm_up_epochs=1,
        save_controller_every=50
    )
    try:
        env.train()
    except KeyboardInterrupt:
        pass
    controller.save_weights(os.path.join(wd, "controller_weights.h5"))
    plot_controller_hidden_states(controller, "%s/controller_states.png" % wd)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Command-line arguments for Search')
    parser.add_argument('--working-dir', "--wd", dest='wd', type=str,
                        help='output working dir')
    parser.add_argument('--GAP', dest='gap', default=False,
                        action="store_true",
                        help='use GlobalAveragePool instead of Flatten')
    parser.add_argument('--layers', dest='layers', type=int, default=12,
                        help='number of layers in model space')
    parser.add_argument('--disable-controller', dest='disable_controller',
                        default=False, action="store_true",
                        help='disable controller training')
    parser.add_argument('--verbose', dest='verbose', default=1,
                        type=int,
                        help='verbose mode')

    args = parser.parse_args()
    if args.wd is not None:
        main(wd=args.wd, num_layers=args.layers, gap=args.gap, disable_controller=args.disable_controller, verbose=args.verbose)
