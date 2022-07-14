library(ggplot2)
setwd("D:/project/Demo_4/SimulationData_Strategy1/")
library(patchwork)
library(R.matlab)
library(SC3)
source("D:/MyWorkWorld/Imputation/scWMC/simulation_1/c_plot_library.r")
# plot the distribution of cell-cell correlation figures in Figure 5 and the supplementary figures realted to Figure 5
# ----------------------------------------------------------------------------
p1 <- list()
for (i in c(1:1)) {
    drop_index <- i
    seed_value <- 10

    p1[[i]] <- plot_cell_distribution(drop_index, seed_value)
}
# -------------------------------------------------------------------------


# plot the distribution of gene-gene correlation figures in Figure 5 and the supplementary figures realted to Figure 5
# ----------------------------------------------------------------------------
p2 <- list()
for (i in c(1:1)) {
    drop_index <- i
    seed_value <- 10

    p2[[i]] <- plot_gene_distribution(drop_index, seed_value)
}

design <- '
ABCDEFGH
IJKLMNOO
'

p = p1[[1]][[1]]
for(i in c(2:15)){
  p = p + p1[[1]][[i]]
}
p = p + plot_layout(guides = 'collect', design = design) + plot_annotation(tag_levels = 'A') & theme(legend.position='bottom')
ggsave(plot = p, filename = "D:/MyWorkWorld/Imputation/scWMC/simulation_1/plots/S1_2.pdf", width = 16, height = 5)

p = p2[[1]][[16]]
for(i in c(17:30)){
  p = p + p2[[1]][[i]]
}
p = p + plot_layout(guides = 'collect', design = design) + plot_annotation(tag_levels = 'A') & theme(legend.position='bottom')
ggsave(plot = p, filename = "D:/MyWorkWorld/Imputation/scWMC/simulation_1/plots/S1_3.pdf", width = 16, height = 5)
