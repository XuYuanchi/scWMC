library(ggplot2)
library(patchwork)
setwd("D:/MyWorkWorld/Imputation/scWMC/ICT")
j=2
data_name = c('Deng', 'Petropoulos')
re = readRDS(paste0("Results/", data_name[j], "/Results.rds"))

methods_name = c("Raw", "scWMC", "ALRA", "DCA", "DrImpute", "EnImpute", "MAGIC", "PBLR", "SAVER", 
                 "scImpute", "scTSSR", "scVI", "VIPER")
p=list()

for(i in c(1:13)){
  p[[i]] = re[[i]][[1]] + 
      ggtitle(methods_name[i],subtitle = paste0("Cor = ", round(abs(re[[i]][[2]]),2), "  POS = ", round(abs(re[[i]][[3]]),2))) 
      # annotate("text", x = -Inf, y = Inf, label = paste0("Cor = ", round(abs(re[[i]][[2]]),2)), hjust = -0.2, vjust = 2, size = 3) + 
      # annotate("text", x = -Inf, y = Inf, label = paste0("POS = ", round(abs(re[[i]][[3]]),2)), hjust = -1.5, vjust = 2, size = 3)

}


design = "
ABCDEFG
HIJKLM#"

pp = p[[1]]
for (i in c(2:13)) {
  pp = pp + p[[i]]
}

pp = pp + plot_layout(guides = "collect", design = design) & 
  theme(plot.title = element_text(size = 8, face = "bold"), plot.subtitle = element_text(size = 8, face = "bold"), legend.position = "bottom")
  
ggsave(plot = pp, filename = paste0("Results/", data_name[j], "/ICT_", data_name[j], ".pdf"), height = 5, width = 14)
