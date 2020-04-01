import sys
import os
from BioNAS.Controller.reward import LossAucReward
from BioNAS.Interpret.scores import PrecisionAtRecall, TprAtFpr
import tensorflow as tf
import keras.backend as K
from src.utils.common_func_msk import read_label_annot, make_output_annot, MotifSaliency, resource_filename
from keras.models import load_model
from src.utils.read_data import read_test_data

s = tf.Session()
K.set_session(s)

label_annot = read_label_annot()
cat_list = [('TF', 'Pol', 'DNase', 'Histone')]
output_annot = make_output_annot(label_annot, cat_list)

def main(model_fp, fn="output_scores.txt"):

    model = load_model(model_fp)
    test_data = read_test_data()

    line_formatter = "{label_name}\t{label_category}\t{fun_name}\t{cutpoint}\t{score}\n"
    with open(fn, "w") as f:
        f.write(line_formatter.replace("{", "").replace("}", "") )

        test_pred = model.predict(test_data[0])
        reward_fns = {
            'PrAtRe': PrecisionAtRecall(cut_points=[0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9], pred=test_pred, compare_op='max'),
            'TpAtFp': TprAtFpr(cut_points=[0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9], pred=test_pred, compare_op='min'),
        }
        for fun_name in reward_fns:
            reward_fn = reward_fns[fun_name]
            auc_list = reward_fn(model, test_data)
            for annot in output_annot:
                scores = auc_list[annot['block']][annot['index']]
                for cutpoint in reward_fn.cut_points:
                    f.write(line_formatter.format(
                                                  fun_name=fun_name, 
                                                  cutpoint=cutpoint,
                                                  score=scores[cutpoint],
                                                  label_name=annot['label_name'],
                                                  label_category=annot['label_category']
                                                  ))
                    f.flush()


if __name__ == '__main__':
    #pass
    model_fp, fn = sys.argv[1], sys.argv[2]
    main(model_fp=model_fp, fn=fn)
