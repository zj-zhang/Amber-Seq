# plotting functions for visualization of Allele-specific DNase and TFs
# ZZ, 2020.2.25

library(ggplot2)
source("darts_theme.R")
library(PRROC)
library(reshape2)
library(dplyr)
library(ggpubr)


plot_scatter = function(data) 
{
	p = ggplot(aes(x=Ref_pred, y=Alt_pred), data=data) + geom_point(aes(colour=ratio-0.5)) +
        scale_colour_gradient2()
	p
}


plot_accuracy_dgf = function(data, margins=seq(0, 0.4, 0.05))
{
	cell_lines = unique(data$new_celltype)
	res_df = NULL
	for(i in 1:length(cell_lines))
	{
		cell_df = NULL
		cell = cell_lines[i]
		this = data[data$new_celltype==cell_lines[i],]
		for(j in 1:length(margins))
		{
			margin = margins[j]
			idx = which( abs(this$pred_margin)>margin )
			if(! length(idx)) next
			acc = mean( this$obs_label[idx]==this$pred_label[idx], na.rm=T )
			cell_df = rbind.data.frame(cell_df, data.frame( task=cell, margin=margin, acc=acc ) )

		}
		res_df = rbind.data.frame(res_df, cell_df)
	}
	
	p = ggplot(aes(x=margin, y=acc, group=task), data=res_df) + geom_line() + ylim(0.5, 1) + Darts_theme

	return(list(df=res_df, fig=p))
}


plot_accuracy_tf = function(data, margins=seq(0, 0.4, 0.05))
{
	tasks = unique(data$deepsea_col)
	res_df = NULL
	for(i in 1:length(tasks))
	{
		task_df = NULL
		task = tasks[i]
		this = data[data$deepsea_col==task,]
		for(j in 1:length(margins))
		{
			margin = margins[j]
			idx = which( abs(this$pred_margin)>margin )
			if(! length(idx)) next
			acc = mean( this$obs_label[idx]==this$pred_label[idx], na.rm=T )
			if(is.nan(acc)) next
			task_df = rbind.data.frame(task_df, data.frame( task=task, margin=margin, acc=acc ) )

		}
		res_df = rbind.data.frame(res_df, task_df)
	}
	
	p = ggplot(aes(x=margin, y=acc, group=task), data=res_df) + geom_line() + ylim(0.5, 1) + Darts_theme

	return(list(df=res_df, fig=p))
}


read_dgf = function(fp)
{
	data = read.table(fp, sep="\t", header=T)
	colnames(data) = c("chr", "pos", "ref", "alt", "ref_reads", "alt_reads", "filename", "celltype", "obs_margin", "ref_pred_2015",
		"alt_pred_2015", "new_celltype", "col_idx", "row_idx", "Ref_pred", "Alt_pred")
	data$total_reads = data$ref_reads + data$alt_reads
	data$ratio = data$ref_reads / data$total_reads
	data$pred_margin = data$Ref_pred - data$Alt_pred
	p_binomial = c()
	for(i in 1:nrow(data)){
		p_binomial = c(p_binomial, binom.test(x=data$ref_reads[i], n=data$total_reads[i], p=0.5)$p.value )
	}
	data$p_binomial = p_binomial
	data$obs_label = NA
	data$obs_label[data$ratio>0.6] = 1
	data$obs_label[data$ratio<0.4] = -1
	data$pred_label = ifelse(data$pred_margin>0, 1, -1)

	data = data[data$p_binomial<0.01, ]
	data = data[data$total_reads>100, ]
	return(data)
}


read_tf = function(fp)
{
	data = read.table(fp, sep="\t", header=T)
	data$ratio = data$ref_reads / data$total_reads
	data$pred_margin = data$Ref_pred - data$Alt_pred
	data$obs_label = NA
	data$obs_label[data$ratio>0.6] = 1
	data$obs_label[data$ratio<0.4] = -1
	data$pred_label = ifelse(data$pred_margin>0, 1, -1)

	data = data[data$p_binomial<0.01, ]
	data = data[data$total_reads>50,]
	return(data)
}


loo_jackknife = function(obs_label, pred_label)
{
	accs = c()
	for(i in 1:length(obs_label))
	{
		this_obs = obs_label[-i]
		this_pred = pred_label[-i]
		accs = c(accs, mean(this_obs==this_pred, na.rm=T))
	}
	return(accs)
}


plot_accuracy_overall = function(data, margins=seq(0,0.4,0.05))
{
	res_df = NULL
	this = data
	for(j in 1:length(margins))
	{
		margin = margins[j]
		idx = which( abs(this$pred_margin)>margin )
		if(! length(idx)) next
		acc = mean( this$obs_label[idx]==this$pred_label[idx], na.rm=T )
		if(is.nan(acc)) next
		
		# for asymptotic ci
		num_total = sum(!is.na(this$obs_label[idx]))
		std = sqrt( acc*(1-acc) / num_total )

		# for jackknife ci
		#accs = loo_jackknife(obs_label=this$obs_label[idx], pred_label=this$pred_label[idx]  )
		#jk_ci = quantile(accs, c(0.025, 0.975))
		#jk_std = sqrt( var(accs) / num_total)

		res_df = rbind.data.frame(res_df, 
					  data.frame( 
						     margin=margin, 
						     acc=acc,
						     num_correct=sum( this$obs_label[idx]==this$pred_label[idx], na.rm=T ),
						     num_total=num_total, 
						     ci_upper=acc+1.96*std,
						     ci_lower=acc-1.96*std
						     #jk_ci_upper=acc+1.96*jk_std,
						     #jk_ci_lower=acc-1.96*jk_std
						     #jk_ci_upper=as.numeric(jk_ci[1]),
						     #jk_ci_lower=as.numeric(jk_ci[2])
						     ) 
		)

	}

	p = ggplot(aes(x=margin, y=acc), data=res_df) + geom_line() + ylim(0.5, 1) + Darts_theme
	return(list(df=res_df, fig=p))
}


plot_comparison = function(df.s, df.r, data_fn="data.tsv")
{
	df.s$category = 'BioNAS'
	df.r$category = 'SampArch'
	
	df = rbind.data.frame(df.s, df.r)
	df$category = factor(df$category, levels=c("search", "random"))
	dodge = 0.015
	p = ggplot(aes(x=margin, y=acc, colour=category), data=df) + geom_point(position=position_dodge(width = dodge), size=1.5) + 
		geom_errorbar(aes(ymin=ci_lower, ymax=ci_upper), width=0, position=position_dodge(dodge)) + 
		scale_colour_manual(values=c("red", "black")) + 
		geom_line(linetype="dashed", position=position_dodge(dodge)) + 
		ylab("Accuracy") + xlab("Prediction Margin") + 
		Darts_theme +
		theme(legend.position=c(1,0), legend.justification=c(1,0), legend.title=element_blank())
	write.table(df, file=data_fn, sep="\t", quote=F, row.names=F )
	return(p)

}


# now use "plot_comparison2" to replace "plot_comparison" to
# restrict the plotting up to 0.1 margin, and add Fisher test p-value
# ZZ, 2020.3.5
plot_comparison2 = function(df.s, df.r, data_fn=NULL)
{

	df.s$category = 'BioNAS'
	df.r$category = 'SampArch'
	
	margin_to_plot = as.factor(c(0, 0.05, 0.1))
	df.s0 = df.s[as.factor(df.s$margin) %in% margin_to_plot,]
	df.r0 = df.r[as.factor(df.r$margin) %in% margin_to_plot,]

	pvals = c()
	for(i in 1:nrow(df.s0))
	{
		cont_table = matrix(c( df.s0$num_correct[i], df.s0$num_total[i], df.r0$num_correct[i], df.r0$num_correct[i] ), nrow=2, ncol=2)
		pvals = c(pvals, fisher.test(cont_table)$p.value)
	}

	df = rbind.data.frame(df.s0, df.r0)
	df$category = factor(df$category, levels=c("BioNAS", "SampArch"))
	#df$margin = factor(df$margin )
	dodge = 0.025
	p = ggplot(aes(x=margin, y=acc, colour=category), data=df) + geom_point(position=position_dodge(width = dodge), size=2) + 
		geom_errorbar(aes(ymin=ci_lower, ymax=ci_upper), width=0.0, position=position_dodge(dodge)) + 
		scale_colour_manual(values=c("red", "black")) + 
		#geom_line(linetype="dashed", position=position_dodge(dodge)) + 
		ylab("Accuracy") + xlab("Prediction Margin") + 
		scale_x_continuous(breaks=c(0, 0.05, 0.1)) + 
		Darts_theme +
		theme(legend.position=c(1,0), legend.justification=c(1,0), legend.title=element_blank()) + 
		annotate("text", x=c(0, 0.05, 0.1), y=1, label=formatC(pvals, format="e", digits=2) )

	if(!is.null(data_fn)) write.table(df, file=data_fn, sep="\t", quote=F, row.names=F )
	return(p)

}


compute_group_auc = function()
{
	fp1 = "./batch_run_figures/3_asb/final_12_5.TF.Allelic_Specific.tsv"
	fp2 = "./batch_run_figures/3_asb/final_12_1_Random.TF.Allelic_Specific.tsv"
	raw = read.table(fp1, sep="\t", header=T)
	raw2 = read.table(fp2, sep="\t", header=T)
	stopifnot( nrow(raw)==nrow(raw2)) 

	data = raw[, c("chr", "pos", "ref", "alt", "tf", "ref_reads", "alt_reads", "total_reads", "fdr_binomial", "deepsea_col")]
	data$ratio = data$ref_reads / data$total_reads
        
	data$nas.search = log(raw$Alt_pred/(1-raw$Alt_pred)) - log(raw$Ref_pred/(1-raw$Ref_pred))
	data$nas.random = log(raw2$Alt_pred/(1-raw2$Alt_pred)) - log(raw2$Ref_pred/(1-raw2$Ref_pred))
	data$deepsea = raw$ds_logfc
	data$deltasvm = raw$deltasvm
	data$deepbind = raw$deepbind_max_delta_raw
	data$jaspar = raw$jaspar_max_delta_raw
	data$meme = raw$best_meme_max_delta_raw


	data$obs_label = NA
        data$obs_label[which(data$ratio>0.6 & data$fdr_binomial<0.01) ] = 1
        data$obs_label[which(data$ratio<0.4 & data$fdr_binomial<0.01) ] = -1
	data$obs_label[data$fdr_binomial>0.9] = 0

        data = data[ which( (data$ref_reads>10 | data$alt_reads>10) & !is.na(data$obs_label)) ,]

	auc_sum_func = function(obs, pred, type="pos"){
		if(type=="pos"){
			scores0 = - pred[which(obs==0 & !is.na(pred) )]
			scores1 = - pred[which(obs==1 & !is.na(pred) )]
		} else {
			scores0 = pred[which(obs==0 & !is.na(pred) )]
			scores1 = pred[which(obs==-1 & !is.na(pred) )]
		}
		if(length(scores1)<20 || length(scores0)<20) {
			return(NA)
		} else {
			auc = roc.curve(scores.class0=scores1, scores.class1=scores0, curve=F)
			return(auc$auc)
		}
	}
	
	tidy_summarize_auc = function(type){
		res = as.data.frame(
		data %>% 
		group_by(deepsea_col) %>%
		summarize(
			  nas.search = auc_sum_func(obs=obs_label, pred=nas.search, type=type),
			  nas.random = auc_sum_func(obs=obs_label, pred=nas.random, type=type),
			  deepsea = auc_sum_func(obs=obs_label, pred=deepsea, type=type),
			  deltasvm = auc_sum_func(obs=obs_label, pred=deltasvm, type=type),
			  deepbind = auc_sum_func(obs=obs_label, pred=deepbind, type=type),
			  jaspar = auc_sum_func(obs=obs_label, pred=jaspar, type=type),
			  meme = auc_sum_func(obs=obs_label, pred=meme, type=type)
		)	
		)
		res0 = as.matrix(res[,-1])
		#print(apply(res0, 2, function(x) mean(x, na.rm=T)) )
		print(apply(res0, 2, function(x) median(x, na.rm=T)) )
		return(res)
	}

	res1 = tidy_summarize_auc(type="pos")
	res2 = tidy_summarize_auc(type="neg")

	res1 = melt(res1)
	p1 = ggplot(aes(x=variable, y=value, colour=variable), data=res1) + 
		geom_boxplot(outlier.shape=NA) + 
		#geom_jitter() + 
		scale_x_discrete(limits = rev(levels(res1$variable))) + 
		scale_y_continuous(breaks=seq(0.4,0.9,0.1), labels=seq(0.4, 0.9, 0.1)) + 
		xlab("") + ylab("AUROC") + 
		coord_flip() + 
		Darts_theme + 
		theme(axis.text.y = element_text(angle=0),  legend.title=element_blank(), legend.position="none") + 
		geom_hline(yintercept=seq(0.4, 0.9, 0.1), linetype="dashed", colour="grey")


	res2 = melt(res2)
	p2 = ggplot(aes(x=variable, y=value, colour=variable), data=res2) + 
		geom_boxplot(outlier.shape=NA) + 
		scale_x_discrete(limits = rev(levels(res2$variable))) + 
		scale_y_continuous(breaks=seq(0.4,0.9,0.1), labels=seq(0.4, 0.9, 0.1)) + 
		xlab("") + ylab("AUROC") + 
		coord_flip() + 
		Darts_theme + 
		theme(axis.text.y = element_text(angle=0),  legend.title=element_blank(), legend.position="none")  + 
		geom_hline(yintercept=seq(0.4, 0.9, 0.1), linetype="dashed", colour="grey")


	ggsave(p1, file="./batch_run_figures/3_asb/TF-gain_comparison.pdf", width=4.5, height=4.5)
	ggsave(p2, file="./batch_run_figures/3_asb/TF-loss_comparison.pdf", width=4.5, height=4.5)

	# updated on 2020.3.5: using `ggpubr`, add p-values
	res1 = tidy_summarize_auc(type="pos")
	res2 = tidy_summarize_auc(type="neg")
	my_comparisons <- list( c("nas.search", "nas.random")
			       #c("nas.search", "deepsea"),
			       #c("nas.search", "deltasvm"),
			       #c("nas.search", "deepbind")
       	)
	idx0.1 = which(apply(is.na(res1[,-1]), 1, sum)!=ncol(res1)-1)
	res0.1 = res1[idx0.1, ]
	res0.1 = melt(res0.1)
	p0.1 = ggboxplot(x="variable", y="value", color="variable", data=res0.1, outlier.shape=NA ) +
		stat_compare_means(comparisons = my_comparisons, method="wilcox.test", paired=T) +
                scale_x_discrete(limits = rev(levels(res0.1$variable)), 
				 labels=rev(c("BioNAS", "SampArch", "DeepSEA", "deltaSVM", "DeepBind", "Jaspar", "MEME"))) +
                scale_y_continuous(breaks=seq(0.4,0.9,0.1), labels=seq(0.4, 0.9, 0.1)) +
                xlab("") + ylab("AUROC") +
                coord_flip() +
                Darts_theme +
                theme(axis.text.y = element_text(angle=0, hjust=1),  legend.title=element_blank(), legend.position="none") +
                geom_hline(yintercept=seq(0.4, 0.9, 0.1), linetype="dashed", colour="grey")

	idx0.2 = which(apply(is.na(res2[,-1]), 1, sum)!=ncol(res2)-1)
	res0.2 = res2[idx0.2, ]
	res0.2 = melt(res0.2)
	p0.2 = ggboxplot(x="variable", y="value", color="variable", data=res0.2, outlier.shape=NA) +
		stat_compare_means(comparisons = my_comparisons, method="wilcox.test", paired=T) +
                scale_x_discrete(limits = rev(levels(res0.1$variable)), 
				 labels=rev(c("BioNAS", "SampArch", "DeepSEA", "deltaSVM", "DeepBind", "Jaspar", "MEME"))) +
                scale_y_continuous(breaks=seq(0.4,0.9,0.1), labels=seq(0.4, 0.9, 0.1)) +
                xlab("") + ylab("AUROC") +
                coord_flip() +
                Darts_theme +
                theme(axis.text.y = element_text(angle=0, hjust=1),  legend.title=element_blank(), legend.position="none") +
                geom_hline(yintercept=seq(0.4, 0.9, 0.1), linetype="dashed", colour="grey")


	ggsave(p0.1, file="./batch_run_figures/3_asb/TF-gain_comparison.withP.pdf", width=4.5, height=4.5)
	ggsave(p0.2, file="./batch_run_figures/3_asb/TF-loss_comparison.withP.pdf", width=4.5, height=4.5)

}



main = function()
{
	data_search.dgf = read_dgf("./batch_run_figures/3_asb/final_12_5.DGF.Allelic_Specific.tsv")
	data_search.tf = read_tf("./batch_run_figures/3_asb/final_12_5.TF.Allelic_Specific.tsv")
	data_random.dgf = read_dgf("./batch_run_figures/3_asb/final_12_1_Random.DGF.Allelic_Specific.tsv")
	data_random.tf = read_tf("./batch_run_figures/3_asb/final_12_1_Random.TF.Allelic_Specific.tsv")

	# task-level individual plots are not stable	
	#res_search.dgf = plot_accuracy_dgf(data_search.dgf)
	#res_search.tf = plot_accuracy_tf(data_search.tf)
	#res_random.dgf = plot_accuracy_dgf(data_random.dgf)
	#res_random.tf = plot_accuracy_tf(data_random.tf)

	# overall accuracy comparison
	res_search.dgf = plot_accuracy_overall(data_search.dgf)
	res_random.dgf = plot_accuracy_overall(data_random.dgf)
	
	res_search.tf = plot_accuracy_overall(data_search.tf)
	res_random.tf = plot_accuracy_overall(data_random.tf)

	#p1 = plot_comparison(res_search.tf$df, res_random.tf$df, data_fn="./batch_run_figures/3_asb/TF.tsv")
	#p2 = plot_comparison(res_search.dgf$df, res_random.dgf$df, data_fn="./batch_run_figures/3_asb/DGF.tsv")
	p1 = plot_comparison2(res_search.tf$df, res_random.tf$df)
	p2 = plot_comparison2(res_search.dgf$df, res_random.dgf$df)
	
	ggsave(p1, file="./batch_run_figures/3_asb/TF_compr2.pdf", width=4.5, height=4.5)
	ggsave(p2, file="./batch_run_figures/3_asb/DGF_compr2.pdf", width=4.5, height=4.5)
}
