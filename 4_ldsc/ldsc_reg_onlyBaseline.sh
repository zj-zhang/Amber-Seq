# shell for running LDSC-reg for Annot
# ZZ, 2020.3.16
# Example call
# $ for i in `cat ./ldsc_resources/ukbb_sumstats_alkesprice/sumstats_list.txt`; do bash ldsc_reg_onlyBaseline.sh $i ; done`


SUMSTAT=$1
DIRNAME=$(basename $SUMSTAT)
mkdir -p ./ldsc_resources/OnlyBaselineResults/$DIRNAME

./ldsc_resources/ldsc/ldsc.py \
	--h2 $SUMSTAT \
	--ref-ld-chr ./ldsc_resources/baselineLD/baselineLD. \
	--frqfile-chr ./ldsc_resources/frq_file/1000G_Phase3_frq/1000G.EUR.QC. \
	--w-ld-chr ./ldsc_resources/weight_file/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
	--overlap-annot \
	--print-coefficients \
	--print-delete-vals \
	--chisq-max 10000 \
	--out ./ldsc_resources/OnlyBaselineResults/$DIRNAME/ldsc_result


