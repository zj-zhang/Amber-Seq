import os
import sys
import pickle
from BioNAS.Controller.model_space import get_layer_shortname
from common_func_msk import read_controller_train_history
import pandas

last_only = int(sys.argv[1]) if len(sys.argv)>1 else None
pattern = "batch_run"
batch_run_dirs = [p for p in os.listdir("./") if os.path.isdir(p) and p.startswith(pattern)]

i = 0
for parDir in batch_run_dirs:
    for subDir in os.listdir(parDir):
        if not subDir.startswith("tmp_search"):
            continue
        model_space_fp = os.path.join(parDir, subDir, "model_space.pkl")
        if not os.path.isfile(model_space_fp):
            continue
        if subDir.endswith("noSearch"):
            is_randomized=True
        else:
            is_randomized=False
        try:
            best_arc, best_auc = read_controller_train_history(os.path.join(parDir, subDir, "train_history.csv"),
                    last_only=last_only,
                    is_randomized=is_randomized)
        except pandas.errors.EmptyDataError:
            continue
        best_arc = ",".join([str(x) for x in best_arc])
        auc_type = "ROC" if best_auc > 0.5 else "PR"
        model_space = pickle.load(open(model_space_fp, "rb"))
        model_states = ",".join([get_layer_shortname(x) for x in model_space[0]]) 
        i += 1
        print("\t".join([str(i), parDir, subDir, auc_type, "%.5f"%best_auc, best_arc, model_states]))
