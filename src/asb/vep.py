from BioNAS.Interpret.sequence_model import AnalyzeSequencesNAS
from selene_sdk.utils import load_features_list
from selene_sdk.sequences import Genome
import argparse
from keras.models import load_model


class NonStrandModel(object):
    def __init__(self, model):
        self.model = model

    def predict_swapbase(self, x):
        x_ = x[:, ::-1, ::-1]
        y_ = self.model.predict(x_)
        return y_

    def predict(self, x):
        y_f = self.model.predict(x)
        y_r = self.predict_swapbase(x)

        return (y_f + y_r)/2.



def analysis(model_path, output_dir, genome_fp, vcf_fp, label_annot_fp):
    label_annot =  load_features_list(label_annot_fp)
    #base_model = load_model(model_path)
    #nsm = NonStrandModel(model=base_model)
    ans = AnalyzeSequencesNAS(
            trained_model_path=model_path,
            #trained_model_path=nsm,
            sequence_length=1000,
            features=label_annot,
            batch_size=1000,
            reference_sequence=Genome(genome_fp),
            swapbase=['A', 'G', 'C', 'T']
            )

    ans.variant_effect_prediction(
            vcf_fp,
            save_data=["abs_diffs", "predictions"],
            output_format="hdf5",
            output_dir=output_dir)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Command-line arguments for Variant Effect prediction')

    parser.add_argument('--model', dest='model', type=str,
                        help='trained keras model path')
    
    parser.add_argument('--od', dest='od', type=str,
                        help='output folder')
    
    parser.add_argument('--genome', type=str,
                        help='genome FASTA filepath')
    
    parser.add_argument('--vcf', type=str,
                        help='VCF filepath')
    
    parser.add_argument('--label-annot', type=str,
                        help='Multi-task label annotation filepath')
    
    parser.add_argument('--verbose', dest='verbose', default=1,
                        type=int,
                        help='verbose mode')

    args = parser.parse_args()
    analysis(
            model_path=args.model, 
            output_dir=args.od,
            genome_fp=args.genome,
            vcf_fp=args.vcf,
            label_annot_fp=args.label_annot
            )
 
