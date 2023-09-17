# calculate correlation 
get.cor.gene <- function(X, Y) {
  sapply(1:nrow(X), function(i) cor(X[i, ], Y[i, ]))
}

get.cor.cell <- function(X, Y) {
  sapply(1:ncol(X), function(i) cor(X[, i], Y[, i]))
}
# plot
plot_comparison <- function(data, ylabels = "Error", color, dropout_rate){
  
  # this function is used to plot the boxplot with the p-values
  # Parameter in the function
  # data : data is a matrix with 6 columns, the first five colums are the errors
  # the last column is the dropout rates
  # ylabels: the y labels shown on the graph
  # ylim_value : the smallest y value where the pvalue is shown
  # h_ylim: the difference between two pvalues shown in the graph. the best ones
  # is the 10%*ylim_vlaue
  
  # extract the data with the first five columns
  dataV0 = data[c(1:11),]
  
  methods = c("drop", "Bis", "ALRA", "DCA", "DeepImpute", "MAGIC", "SAVER", "scImpute", "SCRABBLE", "scScope", "scVI")
  
  dataV1 = data.frame(as.vector(t(dataV0)))
  
  # calculate the dropout rate
  dropout_rate = round(dropout_rate)
  
  # the number of data using for plotting the boxplot
  # dataV1 is built with two columns: y values and group labels
  N = dim(dataV1)[1]     
  
  dataV1$Methods = rep(methods, each = N/11)
  
  colnames(dataV1) = c('y','Methods')
  
  dataV1$Methods = factor(dataV1$Methods, levels=methods)
  
  pp = ggplot(dataV1, aes(x=Methods, y=y, fill=Methods)) +
    geom_boxplot(outlier.shape = NA, position = position_dodge(0.65), width=0.5) +
    #stat_boxplot(geom = "errorbar", width = 0.3)
    scale_fill_manual(values = color) +
    theme_bw() +
    xlab("Methods") + 
    ylab(ylabels) + 
    ggtitle(paste0("Dropout Rate: ",dropout_rate,"%")) +
    theme(axis.title.x=element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(hjust = 0.5),
          text=element_text(size=14),
          axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = rel(0.6)))
  return(pp)
  
}