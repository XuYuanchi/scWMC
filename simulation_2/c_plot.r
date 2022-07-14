library(ggplot2)
setwd("D:/project/Demo_4/SimulationData_Strategy2/")
library(patchwork)
library(R.matlab)
library(SC3)
source("D:/MyWorkWorld/Imputation/scWMC/simulation_2/c_plot_library.r")
# plot the distribution of cell-cell correlation figures in Figure 5 and the supplementary figures realted to Figure 5
# ----------------------------------------------------------------------------
p1 <- list()
for (i in c(1:1)) {
    drop_index <- i
    seed_value <- 1

    p1[[i]] <- plot_cell_distribution(drop_index, seed_value)
}
# -------------------------------------------------------------------------


# plot the distribution of gene-gene correlation figures in Figure 5 and the supplementary figures realted to Figure 5
# ----------------------------------------------------------------------------
p2 <- list()
for (i in c(1:1)) {
    drop_index <- i
    seed_value <- 1

    p2[[i]] <- plot_gene_distribution(drop_index, seed_value)
}

design <- "
ABCDEFG
HIJKLMN
"

p <- p1[[1]][[1]]
for (i in c(2:14)) {
    p <- p + p1[[1]][[i]]
}
p <- p + plot_layout(guides = "collect", design = design) & theme(legend.position = "bottom")
ggsave(plot = p, filename = "D:/MyWorkWorld/Imputation/scWMC/simulation_2/plots/S2_2.pdf", width = 14, height = 5)

p <- p2[[1]][[15]]
for (i in c(16:28)) {
    p <- p + p2[[1]][[i]]
}
p <- p + plot_layout(guides = "collect", design = design) & theme(legend.position = "bottom")
ggsave(plot = p, filename = "D:/MyWorkWorld/Imputation/scWMC/simulation_2/plots/S2_3.pdf", width = 14, height = 5)
