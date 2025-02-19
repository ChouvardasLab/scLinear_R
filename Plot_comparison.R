#Scatterplot of two variables for same observations. I use this to compare predicted ADT values to the base truth
plot_comparison <- function(adt1,adt2,xaxis = 'ADT Values 1', yaxis = 'ADT Values 2',title = ''){
  plot_dat <- cbind(adt1,adt2)
  colnames(plot_dat) <- c('ADT1','ADT2')
  comp_plot <- ggplot(plot_dat,aes(x=ADT1,y=ADT2))+
    geom_point()+
    geom_abline(intercept=0,slope=1,color = "red")+
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    xlab(xaxis) + ylab(yaxis) +
    theme_minimal() +
    coord_equal() +
    ggtitle(title) + theme(
      panel.background = element_rect(fill = "white",
                                      colour = "white",
                                      size = 0.5, linetype = "solid"),
      panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                      colour = "black"),
      panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                      colour = "black")
    )
  return(comp_plot)
}