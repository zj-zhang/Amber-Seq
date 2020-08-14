# plotting functions for visualization of Allele-specific DNase and TFs
# ZZ, 2020.2.25

library(ggplot2)
source("./R/darts_theme.R")
library(PRROC)
library(reshape2)
library(dplyr)
library(ggpubr)
library(grid)
library(gridExtra)
library(ggsci)


# read data and filtering, pre-processing
read_data = function(fp1, fp2)
{
	#fp1 = "./batch_run_figures/3_asb/final_12_5.TF.Allelic_Specific.tsv"
	#fp2 = "./batch_run_figures/3_asb/final_12_1_Random.TF.Allelic_Specific.tsv"
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

	data
}


auc_sum_func = function(obs, pred, type="pos"){
	if(type=="pos"){
		scores0 = - pred[which(obs==0 & !is.na(pred) )]
		scores1 = - pred[which(obs==1 & !is.na(pred) )]
	} else if(type=="neg") {
		scores0 = pred[which(obs==0 & !is.na(pred) )]
		scores1 = pred[which(obs==-1 & !is.na(pred) )]
	}
	else {
		scores0 = abs(pred[which(obs==0 & !is.na(pred) )])
		scores1 = abs(pred[which( abs(obs)==1 & !is.na(pred) )])	
	}
	if(length(scores1)<20 || length(scores0)<20) {
		return(NA)
	} else {
		auc = roc.curve(scores.class0=scores1, scores.class1=scores0, curve=F)
		return(auc$auc)
	}
}


tidy_summarize_auc = function(data, type){
	res = as.data.frame(
	data %>% 
	group_by(deepsea_col) %>%
	summarize(
		  BioNAS = auc_sum_func(obs=obs_label, pred=nas.search, type=type),
		  Sample.Architecture = auc_sum_func(obs=obs_label, pred=nas.random, type=type),
		  DeepSEA = auc_sum_func(obs=obs_label, pred=deepsea, type=type),
		  deltaSVM = auc_sum_func(obs=obs_label, pred=deltasvm, type=type),
		  DeepBind = auc_sum_func(obs=obs_label, pred=deepbind, type=type),
		  Jaspar = auc_sum_func(obs=obs_label, pred=jaspar, type=type),
		  MEME = auc_sum_func(obs=obs_label, pred=meme, type=type)
	),
			    stringsAsFactors=F
	)
	res0 = as.matrix(res[,-1])
	print(apply(res0, 2, function(x) median(x, na.rm=T)) )
	return(res)
}


impute_for_test = function(res)
{
	for(j in 2:ncol(res))
	{
		idx = which(is.na(res[,j]))
		res[idx, j] = mean(res[,j], na.rm=T)
	}
	res
}


# New entry point for plotting, side by side
# 2020.4.21
# ZZ
main = function(fp1, fp2, out)
{
	data = read_data(fp1, fp2)

	# updated on 2020.3.5: using `ggpubr`, add p-values
	
	## LEFT panel: loss-of-func alleles
	res1 = tidy_summarize_auc(data, type="pos")
	idx1 = which(apply(is.na(res1[,-1]), 1, sum)!=ncol(res1)-1)
	res1 = res1[idx1, ]
	res1 = impute_for_test(res1)
	res1.plot = melt(res1)
	res1.plot$variable = factor(res1.plot$variable, levels=rev(c("BioNAS", "Sample.Architecture", "DeepSEA", "deltaSVM", "DeepBind", "Jaspar", "MEME")) )
	g1 = ggboxplot(res1.plot, x="variable", y="value", color="variable", outlier.shape=NA ) +
                coord_flip() +
                #scale_y_reverse(breaks=seq(0.4, 0.9, 0.1), labels=seq(0.4, 0.9, 0.1)) +
                scale_y_continuous(breaks=seq(0.4, 0.9, 0.1), labels=seq(0.4, 0.9, 0.1)) +
		stat_compare_means(label="p.signif", method="wilcox.test", ref.group="BioNAS", paired=T ) +
                xlab("") + ylab("AUROC") +
                Darts_theme +
                theme(axis.text.y = element_text(angle=0, hjust=1),  legend.title=element_blank(), legend.position="none") +
                geom_hline(yintercept=seq(0.4, 0.9, 0.1), linetype="dashed", colour="grey") +
		#scale_color_manual(values=c("BioNAS"="red", "Sample.Architecture"="blue", "DeepSEA"="purple", "deltaSVM"="orange", "DeepBind"="darkgreen", "Jaspar"="darkgrey", "MEME"="brown")) + 
		scale_color_manual(values=c("BioNAS"="red", "Sample.Architecture"="steelblue3", "DeepSEA"="salmon", "deltaSVM"="plum", "DeepBind"="limegreen", "Jaspar"="darkgrey", "MEME"="burlywood")) + 
		ggtitle("Loss-of-function") +
		theme(axis.title.x = element_blank(),
			axis.title.y = element_blank(),
			axis.text.y = element_blank(),
			axis.ticks.y = element_blank(),
			plot.margin = unit(c(1,1.8,1,1), "mm"))


	## RIGHT panel: gain-of-func alleles
	res2 = tidy_summarize_auc(data, type="neg")
	idx2 = which(apply(is.na(res2[,-1]), 1, sum)!=ncol(res2)-1)
	res2 = res2[idx2, ]
	res2 = impute_for_test(res2)
	res2.plot = melt(res2)
	res2.plot$variable = factor(res2.plot$variable, levels=rev(c("BioNAS", "Sample.Architecture", "DeepSEA", "deltaSVM", "DeepBind", "Jaspar", "MEME")) )
	g2 = ggboxplot(res2.plot, x="variable", y="value", color="variable", outlier.shape=NA ) +
		stat_compare_means(label="p.signif", method="wilcox.test", ref.group="BioNAS", paired=T ) +
                coord_flip(ylim=c(0.4, 0.9)) +
                scale_y_continuous(breaks=seq(0.4, 1, 0.1), labels=seq(0.4, 1, 0.1)) +
                xlab("") + ylab("AUROC") +
                Darts_theme +
                theme(axis.text.y = element_text(angle=0, hjust=1),  legend.title=element_blank(), legend.position="none") +
                geom_hline(yintercept=seq(0.4, 0.9, 0.1), linetype="dashed", colour="grey") +
		#scale_color_manual(values=c("BioNAS"="red", "Sample.Architecture"="blue", "DeepSEA"="purple", "deltaSVM"="orange", "DeepBind"="darkgreen", "Jaspar"="darkgrey", "MEME"="brown")) + 
		scale_color_manual(values=c("BioNAS"="red", "Sample.Architecture"="steelblue3", "DeepSEA"="salmon", "deltaSVM"="plum", "DeepBind"="limegreen", "Jaspar"="darkgrey", "MEME"="burlywood")) + 
		ggtitle("Gain-of-function") + 
		theme(axis.title.x = element_blank(),
			axis.title.y = element_blank(),
			axis.text.y = element_blank(),
			axis.ticks.y = element_blank(),
			plot.margin = unit(c(1,2.5,1,0), "mm")) 

	mid.df = res1.plot[which(res1.plot$deepsea_col == res1.plot$deepsea_col[1]),]
	g.mid<-ggplot(mid.df, aes(x=1,y=variable))+geom_text(aes(label=variable))+
	  ggtitle("")+
	  #geom_segment(aes(x=0.94,xend=0.955,yend=variable))+
	  geom_segment(aes(x=1.05,xend=1.065,yend=variable))+
	  ylab(NULL)+
	  theme(axis.title=element_blank(),
		panel.grid=element_blank(),
		axis.text.y=element_blank(),
		axis.ticks.y=element_blank(),
		panel.background=element_blank(),
		axis.text.x=element_text(color=NA),
		axis.ticks.x=element_line(color=NA),
		plot.margin = unit(c(1,-1,1,-1), "mm"))


	#pdf("./batch_run_figures/3_asb/combined_TF.pdf", width=5, height=4.5)
	pdf(out, width=5, height=4.5)
	gg1 <- ggplot_gtable(ggplot_build(g1))
	gg2 <- ggplot_gtable(ggplot_build(g2))
	gg.mid <- ggplot_gtable(ggplot_build(g.mid))
	#grid.arrange(gg1, gg.mid, gg2, ncol=3,widths=c(3/9,2/9,3/9))
	grid.arrange(gg.mid, gg1, gg2, ncol=3,widths=c(2/9,3/9,3/9))
	dev.off()


}


parser = function()
{
        args = commandArgs(T)
        input_search_fp = args[1]
        input_random_fp = args[2]
        output_prefix = args[3]

        # analysis 1
        out = paste0(output_prefix, ".gain_and_loss_.group_auc_withP.pdf")
        main(input_search_fp, input_random_fp, out)

}


parser()

