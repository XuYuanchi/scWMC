plot_comparison <- function(data, ylabels = "Error", ylim_value = 100, h_ylim = 10, color){
  
  # this function is used to plot the boxplot with the p-values
  # Parameter in the function
  # data : data is a matrix with 6 columns, the first five colums are the errors
  # the last column is the dropout rates
  # ylabels: the y labels shown on the graph
  # ylim_value : the smallest y value where the pvalue is shown
  # h_ylim: the difference between two pvalues shown in the graph. the best ones
  # is the 10%*ylim_vlaue
  
  # extract the data with the first five columns
  dataV0 = data[c(1:13),]

  methods = c("Dropout", "scWMC", "ALRA", "DCA", "DrImpute", "EnImpute", "MAGIC", "PBLR", "SAVER",
    "scImpute", "scTSSR", "scVI", "VIPER")
  
  dataV1 = data.frame(as.vector(t(dataV0)))
  
  # calculate the dropout rate
  dropout_rate = round(mean(data[14,]))
  
  # the number of data using for plotting the boxplot
  # dataV1 is built with two columns: y values and group labels
  N = dim(dataV1)[1]     
  
  dataV1$Methods = rep(methods, each = N/13)
  
  colnames(dataV1) = c('y','Methods')
  methods_name = c("Dropout", "scWMC", "ALRA", "DCA", "DrImpute", "EnImpute", "MAGIC", "PBLR", "SAVER", 
                 "scImpute", "scTSSR", "scVI", "VIPER")
  dataV1$Methods = factor(dataV1$Methods, levels=methods_name)

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
  # define the comparison lists
  #my_comparisons <- list( c("1", "13"), c("2", "13"), c("3", "13"), c("4", "13"), c("5", "13"), c("6", "13"), c("7", "13"), c("8", "13"), 
                        #c("9", "13"), c("10", "13"), c("11", "13"), c("12", "13"),)
  
  # compare the errors using t-test
  #pval <- compare_means(y ~ group,data = dataV1, method = "t.test", ref.group = "13", paired = TRUE)
  
  # plot the boxplot with pvalues for the comparisons
#   pp <- ggboxplot(dataV1, x = "group", y = "y", fill = "group",
#                   ) +
#     stat_boxplot(geom = "errorbar", width = 0.3) +
#     #ylim(c(0,ylim_value + 7*h_ylim)) + 
#     theme_bw() +
#     # geom_signif(comparisons = my_comparisons, 
#     #             annotations = formatC(pval$p, format = "e", digits = 2),
#     #             tip_length = 0.01,
#     #             y_position = c(ylim_value + 6*h_ylim, ylim_value + 5*h_ylim, ylim_value + 4*h_ylim, ylim_value + 3*h_ylim, ylim_value + 2*h_ylim, ylim_value + h_ylim, ylim_value)) +
#     theme(text=element_text(size=14)) +
#     xlab("Methods") + 
#     ylab(ylabels) + 
#     ggtitle(paste0("Dropout Rate: ",dropout_rate,"%")) +
#     scale_fill_manual( values = color[1:13]) +
#     theme(panel.grid.major = element_blank(), 
#           panel.grid.minor = element_blank(),
#           panel.background = element_blank())
  return(pp)
  
}