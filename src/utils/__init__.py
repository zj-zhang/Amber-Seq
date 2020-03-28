"""Utility functions for analysis
ZZ 2020.2.24
"""


def deepsea_label_name_normalizer(label_name):
    new_name = label_name.replace("|","--").replace("(","_").replace(")","_")
    return new_name



def run_from_ipython():
    try:
        __IPYTHON__
        return True
    except NameError:
        return False


