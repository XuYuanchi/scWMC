setwd("D:/MyWorkWorld/Imputation/scWMC/ICT")
library(monocle3)
library(gridExtra)
library(ggplot2)
library(scales)
library(TSCAN)
library(R.matlab)

# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds, time_bin="zygote"){
  cell_ids <- which(colData(cds)[, "Time"] == time_bin)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}

process <- function(d, cell_meta, gene.anno, methods_name, cell_label){
cds.raw <- new_cell_data_set(d,
                             cell_metadata = cell_meta,
                             gene_metadata = gene.anno)
#pre-process the data, using PCA with 100 components
cds.raw <- preprocess_cds(cds.raw, num_dim = 5)

# reduction the dimension using one of "tSNE", "PCA", "LSI", "Aligned"
cds.raw <- reduce_dimension(cds.raw, preprocess_method="PCA") # PCA is default, "tSNE", "PCA", "LSI", "Aligned"

cds.raw <- cluster_cells(cds.raw, reduction_method = "UMAP")
cds.raw <- learn_graph(cds.raw, use_partition = F, close_loop = F,
                       learn_graph_control = NULL, verbose = FALSE)

raw.p1 <- plot_cells(cds.raw, color_cells_by="Time", group_label_size = 6, cell_size = 2,
                     label_cell_groups=F,
                     label_leaves=F,
                     label_branch_points=F,
                     graph_label_size=4)+ggtitle(methods_name)+theme(text = element_text(size=8),legend.position = "bottom")

cds.raw <- order_cells(cds.raw, root_pr_nodes=get_earliest_principal_node(cds.raw))

# cor.kendall = cor(cds.raw@phenoData@data$Pseudotime, as.numeric(cds.raw@phenoData@data$timepoint), 
#                   method = "kendall", use = "complete.obs")

cor.kendall = cor(cds.raw@principal_graph_aux@listData$UMAP$pseudotime, as.numeric(cds.raw@colData@listData$Time), 
                       method = "kendall", use = "complete.obs")

lpsorder2 = data.frame(sample_name = cell_meta$Time, State= cds.raw@colData@listData$Time,
                       Pseudotime = cds.raw@principal_graph_aux@listData$UMAP$pseudotime, rank = rank(cds.raw@principal_graph_aux@listData$UMAP$pseudotime))

lpsorder_rank = dplyr::arrange(lpsorder2, rank)

lpsorder_rank$Pseudotime = lpsorder_rank$rank

lpsorder_rank = lpsorder_rank[-4]

lpsorder_rank[1] <- lapply(lpsorder_rank[1], as.character)

subpopulation <- data.frame(cell = cell_meta$Time, sub = as.numeric(cell_label)-1)

POS <- TSCAN::orderscore(subpopulation, lpsorder_rank['sample_name'])[[1]]

Results = list()
Results[[1]] = raw.p1
Results[[2]] = cor.kendall
Results[[3]] = POS
return(Results)
}


# cor.kendall = matrix(nrow = 13, ncol = 1)
# rownames(cor.kendall) = c("Raw", "scWMC", "ALRA", "DCA", "DrImpute", "EnImpute", "MAGIC", "PBLR", "SAVER", 
#                           "scImpute", "scTSSR", "scVI", "VIPER")
# POS = matrix(nrow = 13, ncol = 1)
# rownames(POS) = c("Raw", "scWMC", "ALRA", "DCA", "DrImpute", "EnImpute", "MAGIC", "PBLR", "SAVER", 
#                   "scImpute", "scTSSR", "scVI", "VIPER")

re = list()
# raw
d <- readRDS("data/Deng_raw.rds")
gene_short_name = row.names(d)
cell_label = factor(colnames(d),
                    levels=c('zygote', 'early 2-cell', 'mid 2-cell', 'late 2-cell',
                             '4-cell', '8-cell', '16-cell', 'early blastocyst',
                             'mid blastocyst', 'late blastocyst'))
colnames(d) <- 1:ncol(d)
rownames(d) <- 1:nrow(d)
gene.anno = data.frame(gene_short_name, row.names = rownames(d))
cell_meta = data.frame(cell_label, row.names = colnames(d))
colnames(cell_meta) = "Time"

re[[1]] = process(d, cell_meta, gene.anno, "Raw", cell_label)


#scWMC
d <- readMat("Results/Deng/scWMC.mat")[[1]]
d=t(d)
colnames(d) <- 1:ncol(d)
rownames(d) <- 1:nrow(d)
re[[2]] = process(d, cell_meta, gene.anno, "scWMC", cell_label)

# ALRA
d <- readRDS("Results/Deng/ALRA.rds")
colnames(d) <- 1:ncol(d)
#gene_short_name = row.names(d)
rownames(d) <- 1:nrow(d)
gene.anno = data.frame(gene_short_name, row.names = rownames(d))

re[[3]] = process(d, cell_meta, gene.anno, "ALRA", cell_label)

# DCA
d <- as.matrix(read.table("Results/Deng/DCA.tsv", header = T, sep = "\t", row.names = 1))
colnames(d) <- 1:ncol(d)
rownames(d) <- 1:nrow(d)
gene.anno = data.frame(gene_short_name, row.names = rownames(d))
re[[4]] = process(d, cell_meta, gene.anno, "DCA", cell_label)

# DrImpute
d <- readRDS("Results/Deng/DrImpute.rds")
colnames(d) <- 1:ncol(d)
rownames(d) <- 1:nrow(d)
gene.anno = data.frame(gene_short_name, row.names = rownames(d))
re[[5]] = process(d, cell_meta, gene.anno, "DrImpute", cell_label)

# EnImpute
d <- readRDS("Results/Deng/EnImpute.rds")[[1]]
colnames(d) <- 1:ncol(d)
rownames(d) <- 1:nrow(d)
gene.anno = data.frame(gene_short_name, row.names = rownames(d))
re[[6]] = process(d, cell_meta, gene.anno, "EnImpute", cell_label)

# MAGIC
d <- as.matrix(read.table("Results/Deng/MAGIC.csv", header = T, sep = ",", row.names = 1))
d=t(d)
colnames(d) <- 1:ncol(d)
rownames(d) <- 1:nrow(d)
gene.anno = data.frame(gene_short_name, row.names = rownames(d))
re[[7]] = process(d, cell_meta, gene.anno, "MAGIC", cell_label)

# PBLR
d <- readMat("Results/Deng/PBLR.mat")[[1]]
colnames(d) <- 1:ncol(d)
rownames(d) <- 1:nrow(d)
re[[8]] = process(d, cell_meta, gene.anno, "PBLR", cell_label)

# SAVER
d <- readRDS("Results/Deng/SAVER.rds")
colnames(d) <- 1:ncol(d)
rownames(d) <- 1:nrow(d)
gene.anno = data.frame(gene_short_name, row.names = rownames(d))
re[[9]] = process(d, cell_meta, gene.anno, "SAVER", cell_label)

# scImpute
d <- readRDS("Results/Deng/scImpute.rds")
colnames(d) <- 1:ncol(d)
rownames(d) <- 1:nrow(d)
gene.anno = data.frame(gene_short_name, row.names = rownames(d))
re[[10]] = process(d, cell_meta, gene.anno, "scImpute", cell_label)

# scTSSR
d <- readRDS("Results/Deng/scTSSR.rds")
colnames(d) <- 1:ncol(d)
rownames(d) <- 1:nrow(d)
gene.anno = data.frame(gene_short_name, row.names = rownames(d))
re[[11]] = process(d, cell_meta, gene.anno, "scTSSR", cell_label)

# scVI 
d <- as.matrix(read.table("Results/Deng/scVI.csv", header = F, sep = ","))
d=t(d)
colnames(d) <- 1:ncol(d)
rownames(d) <- 1:nrow(d)
gene.anno = data.frame(gene_short_name, row.names = rownames(d))
re[[12]] = process(d, cell_meta, gene.anno, "scVI", cell_label)

# # VIPER
d <- readRDS("Results/Deng/VIPER.rds")[[1]]
colnames(d) <- 1:ncol(d)
rownames(d) <- 1:nrow(d)
gene.anno = data.frame(gene_short_name, row.names = rownames(d))
re[[13]] = process(d, cell_meta, gene.anno, "VIPER", cell_label)


saveRDS(re, file = "Results/Deng/Results.rds")
