Darts_theme = ggplot2::theme(
		legend.text = ggplot2::element_text(size = 13),
		plot.title = ggplot2::element_text(size=13, face="bold"), 
		axis.title.y = ggplot2::element_text(size=13), 
		axis.title.x = ggplot2::element_text(size=13),
		axis.text.y = ggplot2::element_text(size=13, angle = 90, hjust = 0.5, vjust=0.5), 
		axis.text.x = ggplot2::element_text(size=13, angle=0, hjust=0.5, vjust=0.5),
		legend.background = ggplot2::element_rect(fill = "transparent", colour = "transparent")) +
	ggplot2::theme(
		panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
		panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black")) +
		ggplot2::theme(panel.border = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),
		panel.grid.minor = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"))	


