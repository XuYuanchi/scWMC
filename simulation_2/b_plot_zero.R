library(ggplot2)
library(patchwork)
setwd("D:/MyWorkWorld/Imputation/scWMC/simulation_2/")

source("b_plot_library.r")
color <- readRDS("D:/MyWorkWorld/Imputation/scWMC/color.rds")

error_list <- readRDS(file = "data/error_zero_all.rds")

p1 <- list()

p1[[1]] <- plot_comparison(error_list[[1]], "Log2(Error)", 1400, 140, color)

p1[[2]] <- plot_comparison(error_list[[2]], "Log2(Error)", 1400, 140, color)

p1[[3]] <- plot_comparison(error_list[[3]], "Log2(Error)", 1400, 140, color)

p = (p1[[1]] | p1[[2]] | p1[[3]]) + plot_layout(guides = 'collect') & theme(legend.position='bottom')

ggsave(plot = p, filename = "D:/MyWorkWorld/Imputation/scWMC/simulation_2/zero_error.pdf", height = 4, width = 9)
