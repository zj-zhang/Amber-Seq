# this script reads in the enrichment difference data, and try to
# correlate it to some basic annotations such as mean Chi2
# ZZ
# 2020.3.18

library(dplyr)
library(ggplot2)
library(ggrepel)
source("darts_theme.R")



get_summary_enrichment_median = function(val1, val2, p1, p2)
{
	idx = which( (val1>1 & p1<0.05) | (val2>1 & p2<0.05) )
	res = median( log(val1[idx]) - log(val2[idx]) )
	return(res)
}



read_baseline_ldsc = function()
{
	par_dir = "./ldsc_resources/OnlyBaselineResults/"
	pht_dir = list.dirs(par_dir, recursive=F)
	labels = c()
	baseCoef = c()
	baseCoefZ = c()
	numEnri = c()
	numCoef = c()
	
	for(i in 1:length(pht_dir))
	{
		fp = file.path(pht_dir[i], "ldsc_result.results")
		data = read.table(fp, header=T, sep="\t")
		data$Coef.p = (1 - pnorm(abs(data$Coefficient_z.score)))*2  # two-sided test
		data$Enrichment.q = p.adjust(data$Enrichment_p, method="fdr")
		data$Coef.q = p.adjust(data$Coef.p, method="fdr")

		labels = c(labels, 
			   strsplit(basename(pht_dir[i]), "\\.")[[1]][1] )
		cutoff=0.05

		baseCoefZ = c(baseCoefZ, data$Coefficient_z.score[1])
		baseCoef = c(baseCoef, data$Coefficient[1])
		numEnri = c(numEnri, sum(data$Enrichment.q<cutoff, na.rm=T))
		numCoef = c(numCoef, sum(data$Coef.q<cutoff, na.rm=T))

	}

	df = data.frame(
		name=labels,
		baseCoefZ=baseCoefZ,
		baseCoef=baseCoef,
		numEnri=numEnri,
		numCoef=numCoef,
		stringsAsFactors=F
	)
	return(df)
}


rename_fn = function(x=NULL)
{
	dict = c("body_BMIz"="BMI", 
		"body_HEIGHTz"="Height", 
		"body_WHRadjBMIz"="Waist-Hip ratio",
		"disease_AID_ALL"="AID",
		"disease_ALLERGY_ECZEMA_DIAGNOSED"="Allergy Eczema",
		"disease_ASTHMA_DIAGNOSED"="Asthma",
		"disease_CARDIOVASCULAR"="Cardiovascular",
		"disease_DERMATOLOGY"="Dermatology",
		"disease_HI_CHOL_SELF_REP"="High cholesterol",
		"disease_HYPOTHYROIDISM_SELF_REP"="Hypothyroidism",
		"disease_RESPIRATORY_ENT"="Respiratory & ENT",
		"disease_T2D"="T2D"
	)
	if(is.null(x)) {
		return(dict)
	} else {
		return(dict[as.character(x)])
	}
}


main = function()
{
	res = read.table("./batch_run_figures/4_ldsc/delta_summary.tsv", header=T, sep="\t")

	res_v2.2 = res[which(res$version=="V2.2"),]

	res2 = as.data.frame(
		res_v2.2 %>% group_by(pht) %>% summarise(
			median=get_summary_enrichment_median(
						     val1=model_enrichment_val.1,
						     p1=model_enrichment_P.1,
						     val2=model_enrichment_val.2,
						     p2=model_enrichment_P.2))
		)

	base_info = read.table("./ldsc_resources/ukbb_sumstats_alkesprice/LDSC_basic_info.tsv", header=T, sep="\t")

	colnames(res2) = c("name", "median")
	dm = merge(res2, base_info, by="name")
	dm = dm[!is.na(dm$median),]


	baselineLD_res = read_baseline_ldsc()

	dm = merge(dm, baselineLD_res, by="name")

	no_outlier_df = dm[dm$name!="disease_T2D",]
	lm.res = lm(no_outlier_df$numEnri ~ no_outlier_df$median)
	dm$name = rename_fn(dm$name)
	cor.res = cor.test(no_outlier_df$numEnri, no_outlier_df$median, method='sp')
	cor_txt = paste0("Spearman's rho=", round(cor.res$estimate, 3), ", p=", round(cor.res$p.value, 3), "\n(T2D excluded)")

	p = ggplot(dm, aes(x=median, y=numEnri, label=name)) +
		geom_point() +
		geom_abline(intercept=lm.res$coefficients[1], slope=lm.res$coefficients[2], linetype="dashed", color="red") +
		annotate("label", x=0.7, y=50, label=cor_txt, fill="grey91") +
		geom_text_repel(force=1.5) + 
		xlab("logfc(Enrichment)") + ylab("No. enriched baselineLD annotations") + 
		xlim(-0.5, 1) + 
		#ylim(10, 60) +
		Darts_theme  

	p2 = p + coord_flip() + xlim(-0.45, 1.05) + ylim(10,60)
	
	ggsave(p, file="./batch_run_figures/4_ldsc/delta_enrichment_VS_baselineLD.pdf", width=7, height=7)
	ggsave(p2, file="./batch_run_figures/4_ldsc/delta_enrichment_VS_baselineLD.flip.pdf", width=7, height=7)
	write.table(dm, "./batch_run_figures/4_ldsc/delta_enrichment_VS_baseline.tsv", quote=F, sep="\t", row.names=F)

}

