# for each fraction, plots all bootstrapped Ch1+/Ch2+ ratios (boxplot) and proportion of total concentration (barplot)
plot_ratios_and_proportions = function(fractional.conc, ratio.matrix, expected, grouped, gene.name, output.name) {

  #create a function to make boxplots from 2.5th percentile to 97.5th percentile
  f = function(x) {
    r <- quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
    names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
    r
  }

  # adjust ratio matrix data to be used by ggplot2
  melted.ratio.matrix = melt(ratio.matrix)

  # put proportion of total concentration data for each fraction in the format needed by ggplot2
  num.groups = length(fractional.conc)
  fractional.conc.df = data.frame(rep(1, num.groups), names(fractional.conc), fractional.conc)
  colnames(fractional.conc.df) = c("Var1", "Var2", "value")
  fractional.conc.df$Var2 = factor(rownames(fractional.conc.df), as.character(rownames(fractional.conc.df)))

  #set up the condition parameters for the plot plot
  if (grouped) {
    xlabel = "Grouped Fractions (Number of Ribosomes)"
    tilt = 0
  } else {
    xlabel = NULL
    tilt = 45
  }

  #make the top Ch1+/Ch2+ Ratio boxplot
  plot1 = ggplot(data = melted.ratio.matrix, aes(Var2, value)) +
    stat_summary(fun.data = f, geom="boxplot") +
    geom_abline(intercept = expected, slope = 0, color = "blue") +
    ylab("Ch1+ / Ch2+ Ratio") +
    xlab(xlabel) +
    ggtitle(gene.name) +
    theme(axis.text.x=element_text(angle=tilt, hjust=1), #manually make minimum theme
          plot.title = element_text(size = 14),
          axis.title = element_text(size = 12),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(colour = "black", fill = "white"))

  #make the proportion of total concentration barplot
  plot2 = ggplot(data = fractional.conc.df, aes(Var2, value)) +
    geom_bar(stat = "identity") +
    ylab("Proportion of Total Concentration") +
    theme(axis.text.x=element_text(angle=tilt, hjust=1), #manually make minimum theme
          plot.title = element_text(size = 14),
          axis.title.y = element_text(size = 12),
          axis.title.x = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(colour = "black", fill = "white"))

  grid.newpage()
  pdf(output.name, width = 6, height = 8) #output to pdf
  grid.draw(rbind(ggplotGrob(plot1), ggplotGrob(plot2), size = "last")) #align the two plots vertically
  dev.off()
}
