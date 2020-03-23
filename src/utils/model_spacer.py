"""This script returns constructs a model space,
and should be run prior to running the Snakemake pipeline
ZZ
2020.3.22
"""

from BioNAS.Controller.model_space import State, ModelSpace


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


