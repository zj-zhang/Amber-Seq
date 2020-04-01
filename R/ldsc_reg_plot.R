# read in UKBB gwas LDscore regression results, and plot the some summaries
# ZZJ, 2020.2.6
# revised 2020.2.25: add main for analyzing "./batch_run_20200212-L12-Dilated10"


library(ggplot2)
library(ggrepel)
source("./R/darts_theme.R")

get_phenotype_ldsr = function(par_dir, annot_df=NULL, p.adj.method='bonferroni'){
	sub_dir = list.dirs(par_dir, recursive=F)
	labels = c()
	num_hits_enrichment = c()
	num_hits_coef = c()
	model_enrichment_val = c()
	model_enrichment_P = c()
	model_coef_P = c()
	model_coef_val = c()
	model_Z = c()
	model_enrichment_std = c()
	model_coef_std = c()
	
	# add more infor
	prop_snps = c()
	prop_h2 = c()

	for(d in sub_dir)
	{
		fp = file.path(d, "ldsc_reg.results")
		if(! file.exists(fp)){
			cat("File not found: \n", fp, "\n")
			break
		}
		data = read.table(fp, header=T, sep="\t")
		#data$Coef.p = (1 - pnorm(data$Coefficient_z.score))    # one-sided test: test coef>0
		data$Coef.p = (1 - pnorm(abs(data$Coefficient_z.score)))*2  # two-sided test
		data$Enrichment.q = p.adjust(data$Enrichment_p, method=p.adj.method)
		data$Coef.q = p.adjust(data$Coef.p, method=p.adj.method)

		labels = c(labels, basename(d))
		cutoff = 0.01
		num_hits_enrichment = c(num_hits_enrichment, sum(data$Enrichment.q < cutoff, na.rm=T))
		num_hits_coef = c(num_hits_coef, sum(data$Coef.q < cutoff, na.rm=T))
		model_enrichment_P = c(model_enrichment_P, data[1, "Enrichment.q"])
		model_coef_P = c(model_coef_P, data[1, "Coef.q"])
		model_enrichment_val = c(model_enrichment_val, data[1, "Enrichment"])
		model_Z = c(model_Z, data[1, "Coefficient_z.score"])
		model_enrichment_std = c(model_enrichment_std, data[1,"Enrichment_std_error"])
		model_coef_std = c(model_coef_std, data[1,"Coefficient_std_error"])

		prop_snps = c(prop_snps, data[1, "Prop._SNPs"])
		prop_h2 = c(prop_h2, data[1, "Prop._h2"])

		# read in standardized coefficient value
		tau_star = read.table(file.path(d, "tau_star.txt"), sep="\t", header=T )
		model_coef_val = c(model_coef_val, tau_star$tau_star[1])
	
	}
	df = data.frame(
			label_name=labels, 
			num_hits_enrichment=num_hits_enrichment, 
			num_hits_coef=num_hits_coef,
			model_enrichment_P=model_enrichment_P,
			model_coef_P=model_coef_P,
			model_enrichment_val=model_enrichment_val,
			model_coef_val=model_coef_val,
			model_Z=model_Z,
			model_enrichment_std=model_enrichment_std,
			model_coef_std=model_coef_std,
			prop_snps=prop_snps,
			prop_h2=prop_h2,
		stringsAsFactors=F		
	)
	if(! is.null(annot_df)){
		df = merge(df, annot_df, all.x=T)
	}
	return(df)
}


get_model_ldsr = function(par_dir, out_dir) {
	annot_df = read.table("./ldsc_resources/labels_annot/total_label_idx.csv", header=T, sep="\t")
	pht_dirs = list.dirs(par_dir, recursive=F)
	pht = c()
	enrichment = c()
	ent_annot = c()
	coef = c()
	coef_annot = c()
	for(pd in pht_dirs)
	{
		print(pd)
		d = get_phenotype_ldsr(pd, annot_df=annot_df, p.adj.method='fdr')
		cutoff = 0.01 # if using FDR, lower the cutoff
		enrichment = c(enrichment, sum(d$model_enrichment_P < cutoff, na.rm=T))
		coef = c(coef, sum(d$model_coef_P < cutoff, na.rm=T))
		pht = c(pht, basename(pd))
		ent_annot = c(ent_annot, paste0(d[d$model_enrichment_P<0.05,]$label_name, collapse=","))
		coef_annot = c(coef_annot, paste0(d[d$model_coef_P<0.05,]$label_name, collapse=","))

		plot_pht_volcano(d, fig_title=basename(pd), savefn=file.path(out_dir, paste0(basename(pd),".pdf")))

	}
	df = data.frame(
		pht=pht,
		coef=coef,
		enrichment=enrichment,
		coef_annot=coef_annot,
		ent_annot=ent_annot,
		stringsAsFactors=F
	)
	return(df)
}


plot_pht_volcano = function(df, fig_title="", savefn=NULL)
{
	df$log10_model_coef_P = -log10(df$model_coef_P)
	df$colour = 'grey'
	cutoff = 0.01
	df$colour[df$model_coef_P<cutoff & df$model_coef_val > 0] = "red"
	df$colour[df$model_coef_P<cutoff & df$model_coef_val < 0] = "blue"
	df$show_label = ""
	df$sig_rank = 1000
	df$sig_rank[df$model_coef_val>0] = rank(df$model_coef_P[df$model_coef_val>0])
	df$sig_rank[df$model_coef_val<0] = rank(df$model_coef_P[df$model_coef_val<0])
	idx = which(df$sig_rank < 5 & df$colour!="grey")
	df$show_label[idx] = paste(df$label_idx[idx], df$target[idx], df$cell[idx], sep="|")
	xmax = max(abs(df$model_coef_val))
	p = ggplot(aes(x=model_coef_val, y=log10_model_coef_P, color=colour, label=show_label), data=df) +
		geom_point() + 
		geom_text_repel() + 
		geom_hline(yintercept= - log10(cutoff), linetype="dashed") + 
		geom_vline(xintercept= 0, linetype="dashed") + 
		xlim(-xmax, xmax) +
		xlab("Standardized effect size") + 
		ylab("-log10(FDR)") + 
		ggtitle(fig_title) + 
		scale_colour_manual(values=c("red"="red", "blue"="blue", "grey"="grey")) + 
		Darts_theme
	
	if(!is.null(savefn)) {
		ggsave(p, file=savefn)
		write.table(df, file.path(dirname(savefn), paste0(basename(savefn), ".data")), quote=F, row.names=F)
	}
	return(p)
}


plot_model_wise_hits = function(res_df, annot_df, savefn=NULL) {
	plot_df = NULL
	for(i in 1:nrow(res_df)){
		this = strsplit(res_df[i,]$coef_annot, ",")[[1]]
		annot = annot_df[annot_df$label_name %in% this, ]
		cat_freq = cbind( "GWAS"=res_df$pht[i], as.data.frame(table(annot$category)))
		colnames(cat_freq) = c("GWAS", "Category", "Count")
		plot_df = rbind.data.frame(plot_df, cat_freq)
	}

	p = ggplot(aes(x=GWAS, y=Count, fill=Category), data=plot_df) + 
		geom_bar(stat="identity", width=0.9, colour="black") + 
		scale_fill_brewer(palette = 'Set1') + 
		Darts_theme +
		theme(plot.title=element_text(hjust=0.5), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1),
			axis.text.x = element_text(size=10, angle=45, hjust=1, vjust=1),axis.text.y = element_text(size=10))
	if(! is.null(savefn))
	{
		ggsave(p, file=savefn)
	}
	return(p)
}


compare_baseline = function(par_dir1, par_dir2, pht)
{
	read_data = function(par_dir) {
		baselines = c("V2.2"="label_wise_ldsc_reg", "V1.1"="label_wise_ldsc_reg_BaselineV1.1", "OnlyLDsize"="label_wise_ldsc_reg_noBaseline_withLDblock")
		res = NULL
		for(i in 1:length(baselines))
		{
			pht_dir = file.path(par_dir, baselines[i], pht)
			pht_df = get_phenotype_ldsr(pht_dir, p.adj.method="fdr" )
			pht_df = pht_df[, c("label_name", "model_enrichment_P", "model_enrichment_val", "model_coef_P", "model_coef_val", "model_enrichment_std", "model_coef_std")]
			pht_df$version = names(baselines)[i]
			res = rbind.data.frame(res, pht_df)
		}

		return(res)
	}
	
	res1 = read_data(par_dir1)
	res2 = read_data(par_dir2)
	
	delta = merge(res1, res2, by=c("label_name", "version"), suffix=c(".1", ".2"))
	delta$model_enrichment_val = delta$model_enrichment_val.1 - delta$model_enrichment_val.2
	delta$model_coef_val = delta$model_coef_val.1 - delta$model_coef_val.2

	delta$pht = pht
	return(delta)
}


main = function()
#compare_multiple_phenotypes = function()
# the new main entry; 2020.3.11
{
	res = NULL
	par_dir1 = "./batch_run_20200212-L12-Dilated10/tmp_final_12_5/vep_output/"
	par_dir2 = "./batch_run_20200212-L12-Dilated10/tmp_final_12_1.Random/vep_output/"
	
	phenotypes = c("body_BMIz", "body_HEIGHTz", "body_WHRadjBMIz", "disease_AID_ALL", "disease_ALLERGY_ECZEMA_DIAGNOSED", "disease_ASTHMA_DIAGNOSED", "disease_CARDIOVASCULAR", "disease_DERMATOLOGY", "disease_HI_CHOL_SELF_REP", "disease_HYPOTHYROIDISM_SELF_REP", "disease_RESPIRATORY_ENT", "disease_T2D")

	if(file.exists("./batch_run_figures/4_ldsc/delta_summary.tsv"))
	{
		res = read.table("./batch_run_figures/4_ldsc/delta_summary.tsv", header=T, sep="\t")
	} else {
		for(i in 1:length(phenotypes))
		{
			pht = phenotypes[i]
			print(pht)
			delta = compare_baseline(par_dir1, par_dir2, pht)
			res = rbind.data.frame(res, delta)
		}

		write.table(res, file="./batch_run_figures/4_ldsc/delta_summary.tsv", quote=F, sep="\t", row.names=F)
	}


	# rename function for plotting
	#rename_fn =function(x){ a=strsplit(as.character(x), "_",)[[1]]; a=a[2:length(a)]; a=paste(a, collapse="_"); return(a)}

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

	# PLOT enrichment
	res1 = res[which(
			 (res$model_enrichment_val.1>1 & res$model_enrichment_P.1<0.05) | 
			 (res$model_enrichment_val.2>1 & res$model_enrichment_P.2<0.05)
			 ), ]
	res1 = res1[ which(
			   res1$version %in% c( "V2.2")  
		   ),]

	res1$model_enrichment_val = log(res1$model_enrichment_val.1) - log(res1$model_enrichment_val.2)
	res1$category = "Disease"
	res1$category[which(res1$pht %in% c("body_BMIz", "body_HEIGHTz", "body_WHRadjBMIz"))] = "Trait"	

	res1$pht2 = rename_fn(res1$pht)
	summary_data = read.table("./batch_run_figures/4_ldsc/delta_enrichment_VS_baseline.tsv", header=T, sep="\t")
	res1$pht2 = factor(res1$pht2, levels = rename_fn() )
	p1 = ggplot(res1, aes(x=pht2, y=model_enrichment_val, fill=category)) + 
		geom_boxplot(outlier.shape=NA, position=position_dodge(0.8)) +
		geom_hline(yintercept=0, linetype="dashed", colour="grey") + 
		#geom_jitter() + 
		ylim(-1, 1.5) + 
		Darts_theme +
		theme( axis.text.x = element_text(size=10, angle=45, hjust=1, vjust=1)) +
	        ylab("logfc(Enrichment)") + xlab("") + 
		theme(legend.position=c(0,1), legend.justification=c(0,1), legend.title=element_blank())
	
	# DO NOT RUN: this will mix across categories, and because negative ones 
	# have more in numbers, the distributions look bad
	# update 2020.3.19: version2, plot histogram
	#p1.2 = ggplot(res1, aes(x=model_enrichment_val, fill=category)) +
	#	      geom_histogram(bins=100)
	#ggsave(p1.2, file="./batch_run_figures/4_ldsc/normal_disease.delta.pdf", width=4.5, height=4.5)


	# PLOT coef
	res2 = res[which(
			 (res$model_coef_val.1>0 & res$model_coef_P.1<0.05) | 
			 (res$model_coef_val.2>0 & res$model_coef_P.2<0.05 )
		 ), ]
	res2 = res2[ which(
			   res2$version %in% c("V2.2")
		   ),]


	res2$category = "Disease"
	res2$category[which(res2$pht %in% c("body_BMIz", "body_HEIGHTz", "body_WHRadjBMIz"))] = "Trait"	
	res2$pht2 = rename_fn(res2$pht)
	res2$pht2 = factor(res2$pht2, levels=rename_fn())

	p2 = ggplot(res2, aes(x=pht2, y=model_coef_val, fill=category)) + 
		geom_boxplot(outlier.shape=NA, position=position_dodge(0.8)) +
		geom_hline(yintercept=0, linetype="dashed", colour="grey") + 
		Darts_theme +
		theme( axis.text.x = element_text(size=10, angle=45, hjust=1, vjust=1)) + 
		ylim(-0.3, 0.4) +
		ylab("Coefficient Margin") + xlab("") + 
		theme(legend.position=c(0,1), legend.justification=c(0,1), legend.title=element_blank())

	# SAVE
	ggsave(p1, file="./batch_run_figures/4_ldsc/delta_enrichment.pdf", width=4.5, height=7)
	ggsave(p2, file="./batch_run_figures/4_ldsc/delta_coef.pdf", width=4.5, height=7)
	write.table(res1, "./batch_run_figures/4_ldsc/delta_enrichment.tsv", sep="\t", quote=F, row.names=F)
	write.table(res2, "./batch_run_figures/4_ldsc/delta_coef.tsv", sep="\t", quote=F, row.names=F)


}


