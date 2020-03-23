import tensorflow as tf
import numpy as np
import keras.backend as K
import os
from BioNAS.Controller.model_builder import EnasCnnModelBuilder
from BioNAS.Controller.general_controller import GeneralController
from BioNAS.Controller.model_space import State, ModelSpace
from BioNAS.Controller.common_ops import count_model_params

from src.utils.read_data import read_val_data, read_train_data, read_test_data, read_label_annot
from keras.callbacks import ModelCheckpoint, EarlyStopping
from keras.optimizers import Adam, SGD

from BioNAS.Controller.manager import EnasManager
from BioNAS.Controller.train_env import EnasTrainEnv, setup_logger
from BioNAS.Controller.reward import LossAucReward
from BioNAS.utils.plots import plot_controller_hidden_states
import logging
import pickle
import json


def get_controller(controller_config, model_space, session):
    with tf.device("/cpu:0"):
        lr = controller_config['learning_rate']
        if not controller_config['with_input_blocks']:
            assert  controller_config['num_input_blocks'] == 1
        controller = GeneralController(
            model_space,
            session=session,
            #share_embedding={i:0 for i in range(1, len(model_space))},
            share_embedding={int(k): int(v) for k,v in controller_config['share_embedding'].items()},
            with_skip_connection=controller_config['with_skip_connection'],
            with_input_blocks=controller_config['with_input_blocks'],
            num_input_blocks=controller_config['num_input_blocks'],
            skip_connection_unique_connection=controller_config['skip_connection_unique_connection'],
            skip_weight=controller_config['skip_weight'],
            skip_target=controller_config['skip_target'],
            lstm_size=controller_config['lstm_size'],
            lstm_num_layers=controller_config['lstm_num_layers'],
            kl_threshold=controller_config['kl_threshold'],
            train_pi_iter=controller_config['train_pi_iter'],
            optim_algo='adam',
            temperature=controller_config['temperature'],
            lr_init=lr,
            lr_dec_rate=1.0,
            tanh_constant=controller_config['tanh_constant'],
            buffer_size=controller_config['buffer_size'],
            batch_size=controller_config['batch_size']
        )
        controller.buffer.rescale_advantage_by_reward = False
    return controller



def get_controller_config(controller):
    d = controller.__dict__
    d = {k:d[k] for k in d  if type(d[k]) not in (tf.Tensor, tf.Variable, list, dict)}
    non_picklable_keys = ('optimizer', 'train_op', 'session', 'g_emb', 'w_attn_1', 'w_attn_2', 'v_attn', 'train_step', 'optim_algo' )
    d = {k:d[k] for k in d if not k in non_picklable_keys}
    print(d)
    return d


def main(wd, model_space, controller_config, gap=False, verbose=0):
    config = tf.ConfigProto()
    config.gpu_options.allow_growth = True
    session = tf.Session(config=config)
    K.set_session(session)
    out_filters = controller_config['child_base_filters']
    controller = get_controller(controller_config, model_space, session=session)
    pickle.dump(model_space, open(os.path.join(wd, "model_space.pkl"), "wb"))
    controller_config_dict = get_controller_config(controller)
    pickle.dump(controller_config_dict, open(os.path.join(wd, "controller_config.pkl"), "wb") )

    #reward_fn = LossAucReward(method='aupr')
    reward_fn = LossAucReward(method='auc')

    input_node = State('input', shape=(1000, 4), name="input", dtype=tf.float32)
    output_node = State('dense', units=919, activation='sigmoid')
    model_compile_dict = {
        'loss': 'binary_crossentropy',
        'optimizer': 'adam',
        #'optimizer': SGD(lr=0.1, momentum=0.9, decay=1e-6, nesterov=False),
    }

    child_batch_size = controller_config['child_batch_size']
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
        max_episode=controller_config['manager_max_episode'],
        max_step_per_ep=controller_config['manager_max_step_per_ep'],
        working_dir=wd,
        time_budget=controller_config['manager_time_budget'],
        with_input_blocks=controller_config['with_input_blocks'],
        with_skip_connection=controller_config['with_skip_connection'],
        child_train_steps=controller_config['child_train_steps'],
        child_warm_up_epochs=controller_config['child_warm_up_epochs'],
        save_controller_every=controller_config['save_controller_every']
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
    
    parser.add_argument('--model-space', dest="model_space", type=str,
                        help='filepath to BioNAS model space')
    
    parser.add_argument('--controller-config', dest="controller_config", type=str,
                        help='Conriguration file for controller')
    
    parser.add_argument('--GAP', dest='gap', default=False,
                        action="store_true",
                        help='use GlobalAveragePool instead of Flatten')
    
    parser.add_argument('--verbose', dest='verbose', default=1,
                        type=int,
                        help='verbose mode')

    args = parser.parse_args()
    # load data files
    assert os.path.isfile(args.controller_config)
    controller_config = json.load(open(args.controller_config, "r"))

    assert os.path.isfile(args.model_space)
    model_space = pickle.load(open(args.model_space, "rb"))

    if args.wd is not None:
        main(wd=args.wd, model_space=model_space, controller_config=controller_config, gap=args.gap, verbose=args.verbose)
