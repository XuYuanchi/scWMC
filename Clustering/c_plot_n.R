# plot ARI and NMI
library(ggplot2)
library(datasets)
library(dplyr)
library(tidyr)
library(patchwork)

methods_name = c("Raw", "scWMC", "ALRA", "DCA", "DrImpute", "EnImpute", "MAGIC", "PBLR", "SAVER", 
                 "scImpute", "scTSSR", "scVI", "VIPER")
dataset_name = c("sc_CELseq2", "sc_10x", "sc_Droseq", "sc_10x_5cl", "sc_Celseq2_5cl_p1", "sc_Celseq2_5cl_p2", "sc_Celseq2_5cl_p3")
color = readRDS("D:/MyWorkWorld/Imputation/scWMC/color.rds")

read_NMI <- function(i, j, methods_name, dataset_name, clustering_name){
  dir = paste0("D:/MyWorkWorld/Imputation/scWMC/Clustering/Re/", clustering_name, 
               "/NMI_")
  file_name = paste0(dir, i, ".rds")
  data = readRDS(file_name)
  nmi = data[, j]
  NMI = data.frame("methods" = methods_name,
                   "dataname" = rep(dataset_name, time=13),
                   "num" = rep(i,time=13))
  NMI$nmi = nmi
  return (NMI)
}

clustering_name = c("Seurat", "SC3")
p1 <- vector("list", 2)
for(m in c(1:1)){
  for(j in c(3:6)){
    NMI = data.frame()
    for(i in c(1:10)){
      NMI = rbind(NMI, read_NMI(i, j, methods_name = methods_name, 
                                dataset_name = dataset_name[j], clustering_name=clustering_name[m]))
    }
    NMI_sum <- group_by(NMI, methods) %>%
      summarise(
        mean = mean(nmi), 
        se = sd(nmi) / sqrt(n()),
        lbar = mean - sd(nmi),
        hbar = mean + sd(nmi))
    
    NMI_sum$methods <- factor(NMI_sum$methods, levels=methods_name)
    
    p1[[m]][[j]] <- ggplot(data = NMI_sum, mapping = aes(x = methods, y = mean, fill=methods)) + 
      geom_bar(stat = "identity") +
      scale_fill_manual(values = color) +
      geom_errorbar(mapping = aes(
        ymin = lbar, ymax = hbar), 
        position = position_dodge(width = 1),
        width = 0.1, size = 0.5) +
      labs(title=dataset_name[j], y = "NMI") +
      #coord_cartesian(ylim=c(0,1)) +
      theme_bw() +
      theme(
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5))
  }
  
}
saveRDS(p1, file = "D:/MyWorkWorld/Imputation/scWMC/Clustering/Re/NMI_Seurat.rds")

p11 = (p1[[1]][[3]] | p1[[1]][[4]] | p1[[1]][[5]] | p1[[1]][[6]]) + 
  plot_layout(guides = "collect") + plot_annotation(tag_levels = 'A') & theme(legend.position = "right", legend.key.size = unit(15,"pt"))

