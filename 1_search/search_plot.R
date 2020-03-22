# read in search data and metrices generated for an individual run, and plot figures
# ZZ
# 2020.2.27

library(ggplot2)
library(ggrepel)
source("darts_theme.R")

plot_prob_over_layers = function()
{
	data = read.table("./batch_run_figures/1_search/selection_prob_at_convergence.tsv", header=T, sep="\t")
	data$layer_type = factor(data$layer_type, levels=c("conv_k8_relu", "conv_k4_relu", "conv_k8_relu_d10", "conv_k4_relu_d10", "maxpool1d", "avgpool1d", "identity"))
	data$layer_id = data$layer_id + 1
	p = ggplot(aes(x=layer_id, y=prob, colour=layer_type), data=data) + 	
		#annotate("rect", xmin=c(0.5, 6.5), xmax=c(3.5, 9.5), ymin=0, ymax=1, alpha=0.2) + 
		annotate("rect", xmin=c(3.5, 9.5), xmax=c(6.5, 12.5), ymin=0, ymax=1, alpha=0.1) + 
		geom_point(size=2) +
		geom_line(linetype="dashed") + 
		xlab("Layers") + ylab("Prob at Convergence") + 
	        xlim(0, 13) + ylim(0,1) +	
		Darts_theme +
		theme(legend.position=c(0,1), legend.justification=c(0,1), legend.title = element_blank())
	ggsave(p, file="./batch_run_figures/1_search/selection_prob_at_conv.pdf", width=4.5, height=4.5)
	return(p)
}


plot_pca = function()
{
	data = read.table("./batch_run_figures/1_search/pca_1-2_data.txt", header=T, sep="\t")

	p = ggplot(aes(x=PC1..0.59, y=PC2..0.17, label=labels, colour=layer_types), data=data) + 
		geom_point(size=2) +
		#geom_label() + 
		geom_text_repel() +
		geom_hline(yintercept=0, linetype="dashed") +
		geom_vline(xintercept=0, linetype="dashed") +
		Darts_theme + 
		theme(legend.position=c(1,1), legend.justification=c(1,1), legend.title = element_blank())
	ggsave(p, file="./batch_run_figures/1_search/PCA_1-2_embed.pdf", width=4.5, height=4.5)
	
	data = read.table("./batch_run_figures/1_search/pca_2-3_data.txt", header=T, sep="\t")

	p = ggplot(aes(x=PC2..0.17, y=PC3..0.11, label=labels, colour=layer_types), data=data) + 
		geom_point(size=2) +
		#geom_label() + 
		geom_text_repel() +
		geom_hline(yintercept=0, linetype="dashed") +
		geom_vline(xintercept=0, linetype="dashed") +
		Darts_theme +
		theme(legend.position=c(1,1), legend.justification=c(1,1), legend.title = element_blank())
	ggsave(p, file="./batch_run_figures/1_search/PCA_2-3_embed.pdf", width=4.5, height=4.5)

}
