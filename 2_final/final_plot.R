# plotting beatified figures for 1_final in R
# ZZ 2020.02.22

library(ggplot2)
source("darts_theme.R")

plot_summary_performance = function()
{
	od = "./batch_run_figures/2_final"
	data = read.table(file.path(od, "raw.tsv"), header=T, sep="\t")


	# plot Precision-Recall
	p = ggplot(data, aes(x=rank, y=tests_aupr, group=category, color=category)) + geom_point(aes(size=params)) + 
		geom_line(aes(group=category), linetype="dashed") + 
		#geom_hline(yintercept=0.338, linetype="solid", colour="black") + 
		scale_colour_manual(values=c("orange", "darkblue"), labels=c("SampArch", "BioNAS")) + 
		#annotate("text", x=1, y=0.338 + 0.002, label="2015-DeepSEA", hjust=0) + 
		xlab("Runs") + ylab("Test AUPR") + ylim(0.3, 0.375) +
		Darts_theme #+ 
		#theme(legend.position=c(1,0), legend.justification=c(1,0))

	ggsave(p, file=file.path(od, "comparison_aupr.pdf"), width=5, height=4.5)
	
	# plot AUROC
	p0 = ggplot(data, aes(x=rank, y=tests_auc, group=category, color=category)) + geom_point(aes(size=params)) + 
		geom_line(aes(group=category), linetype="dashed") + 
		#geom_hline(yintercept=0.931, linetype="solid", colour="black") + 
		scale_colour_manual(values=c("orange", "darkblue"), labels=c("SampArch", "BioNAS")) + 
		#annotate("text", x=2, y=0.931 + 0.0005, label="2015-DeepSEA", hjust=0) + 
		xlab("Runs") + ylab("Test AUROC") + ylim(0.915, 0.94) +
		Darts_theme #+ 
		#theme(legend.position=c(1,0), legend.justification=c(1,0))

	ggsave(p0, file=file.path(od, "comparison_auroc.pdf"), width=5, height=4.5)

}


plot_delta_performance = function()
{
	d1 = read.table("./batch_run_figures/2_final/scorer_output/tmp_final_12_5.scorer_output.txt", header=T)
	d0 = read.table("./batch_run_figures/2_final/scorer_output/tmp_final_12_1.Random.scorer_output.txt", header=T)
	# plot Precision at Recall
	d.delta = d1
	d.delta$score = d1$score - d0$score

	d.delta = d.delta[which(d.delta$fun_name=="PrAtRe" & !is.na(d.delta$score)) ,]
	d.delta$cutpoint = factor(d.delta$cutpoint, 
				  levels=c("0.01", "0.05", "0.1", "0.2", "0.3", "0.4", "0.5",
					   "0.6", "0.7", "0.8", "0.9"
					   ))
	p = ggplot(aes(x=cutpoint, y=score), data=d.delta) + 
		geom_boxplot(outlier.shape=NA) + 
		stat_summary(fun.y=mean, geom="point", shape=23, size=4) +
		ylim(-0.1, 0.15) +
		geom_hline(yintercept=0, linetype="dashed", colour="grey") +
		xlab("Recall") + ylab("Precision difference") +
		Darts_theme

	ggsave(p, file="./batch_run_figures/2_final/PrRe_delta.pdf", width=4.5, height=4.5)

	# plot TPR at FPR
	d.delta2 = d1
	d.delta2$score = d1$score - d0$score

	d.delta2 = d.delta2[which(d.delta2$fun_name=="TpAtFp" & !is.na(d.delta2$score)) ,]
	d.delta2$cutpoint = factor(d.delta2$cutpoint, 
				  levels=c("0.01", "0.05", "0.1", "0.2", "0.3", "0.4", "0.5",
					   "0.6", "0.7", "0.8", "0.9"
					   ))
	p2 = ggplot(aes(x=cutpoint, y=score), data=d.delta2) + 
		geom_boxplot(outlier.shape=NA) + 
		stat_summary(fun.y=mean, geom="point", shape=23, size=4) +
		ylim(-0.05, 0.1) +
		geom_hline(yintercept=0, linetype="dashed", colour="grey") + 
		xlab("FPR") + ylab("TPR difference") +
		Darts_theme

	ggsave(p2, file="./batch_run_figures/2_final/TprFpr_delta.pdf", width=4.5, height=4.5)
}
