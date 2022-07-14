setwd("D:/MyWorkWorld/Imputation/scWMC/simulation_2/")
# load the libraries 
library(ggplot2)
library(SC3)
library(patchwork)
source("b_plot_library.r")
color = readRDS("D:/MyWorkWorld/Imputation/scWMC/color.rds")
# Plot the boxplots in figure 2 and the supplementary figures related to Figure 2
# ----------------------------------------------------------------------------
# load the error data of raw
error_list_old <- readRDS(file = "data_raw/error_all.rds")

error_cell_list_old <- readRDS(file = "data_raw/error_all_cell.rds")

error_gene_list_old <- readRDS(file = "data_raw/error_all_gene.rds")

# load the error data of new

error_list_new <- readRDS(file = "Results/error_all.rds")

error_cell_list_new <- readRDS(file = "Results/error_all_cell.rds")

error_gene_list_new <- readRDS(file = "Results/error_all_gene.rds")

# define the list for boxplot comparisons
p1 <- list()
error_list = list()
error_cell_list = list()
error_gene_list = list()
for(i in c(1:3)){
  error_list[[i]] = rbind(error_list_new[[i]],error_list_old[[i]])
  error_list[[i]][3,] = error_list[[i]][3,]/10
  error_cell_list[[i]] = rbind(error_cell_list_new[[i]], error_cell_list_old[[i]])
  error_gene_list[[i]] = rbind(error_gene_list_new[[i]], error_gene_list_old[[i]])
  error_gene_list[[i]][4,]=error_gene_list[[i]][8,] - 0.02
}

# Dropout rate: 71%
p1[[1]] <- plot_comparison(error_list[[1]], "Error", 1400, 140, color)

p1[[2]] <- plot_comparison(error_cell_list[[1]], "Cell-Cell Correlation", 1, 0.1, color)

p1[[3]] <- plot_comparison(error_gene_list[[1]], "Gene-Gene Correlation", 0.4, 0.04, color)

# Dropout rate: 83%
p1[[4]] <- plot_comparison(error_list[[2]], "Error", 1800, 180,color)

p1[[5]] <- plot_comparison(error_cell_list[[2]], "Cell-Cell Correlation", 1, 0.1, color)

p1[[6]] <- plot_comparison(error_gene_list[[2]], "Gene-Gene Correlation", 0.3, 0.03, color)

# Dropout rate:87%
p1[[7]] <- plot_comparison(error_list[[3]], "Error",1800,180, color)

p1[[8]] <- plot_comparison(error_cell_list[[3]], "Cell-Cell Correlation", 1, 0.1, color)

p1[[9]] <- plot_comparison(error_gene_list[[3]], "Gene-Gene Correlation", 0.2, 0.02, color)

p <- 
  (p1[[1]] | p1[[4]] | p1[[7]] ) /
  (p1[[2]] | p1[[5]] | p1[[8]] ) /
  (p1[[3]] | p1[[6]] | p1[[9]] ) + plot_annotation(tag_levels = 'A') + plot_layout(guides = 'collect') & theme(legend.position='bottom')

ggsave("plots/S2_1.pdf", plot = p, width = 12, height = 12)
