# clustering 
library(Seurat)
library(aricode)
library(stats)
library(R.matlab)
library(umap)
dirc = "D:/MyWorkWorld/Imputation/scWMC/scalability/Results/"
me = list()

cal_metric <- function(count){
  
  ncluster = length(unique(colnames(count)))
  
  km.label.k <- kmeans(t(count), centers = ncluster, nstart = 50)$cluster
  
  label <- colnames(count)
  
  #ARI <- adjustedRandIndex(as.factor(label), x.seurat$seurat_clusters)
  ari <- ARI(as.factor(label), km.label.k)
  nmi <- NMI(as.factor(label), km.label.k)
  metric <- c(ari, nmi)
  return (metric)
}
cal_metric_Seurat <- function(count){
  x.seurat <- CreateSeuratObject(count)
  x.seurat <- NormalizeData(x.seurat)
  x.seurat <- FindVariableFeatures(x.seurat, verbose = FALSE)
  x.seurat <- ScaleData(x.seurat)
  
  x.seurat <- RunPCA(x.seurat, features = VariableFeatures(object = x.seurat))  
  x.seurat <- JackStraw(x.seurat, num.replicate = 100)
  x.seurat <- ScoreJackStraw(x.seurat, dims = 1:20)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
  x.seurat <- FindNeighbors(x.seurat, dims = 1:10)
  x.seurat <- FindClusters(x.seurat, resolution = 0.1)
  
  label <- sub('.*:','',colnames(count))
  
  #ARI <- adjustedRandIndex(as.factor(label), x.seurat$seurat_clusters)
  ari <- ARI(as.factor(label), x.seurat$seurat_clusters)
  nmi <- NMI(as.factor(label), x.seurat$seurat_clusters)
  metric <- c(ari, nmi)
  return (metric)
}
#raw
dat = as.matrix(read.table(paste0("D:/MyWorkWorld/Imputation/scWMC/scalability/data/100k.csv"), sep = ",", header = F, row.names = 1))
dat = readRDS("D:/MyWorkWorld/Imputation/scWMC/scalability/data/100k.rds")
# cell_name = colnames(dat)
# saveRDS(cell_name, file = "D:/MyWorkWorld/Imputation/scWMC/scalability/data/cell_name.rds")

cell = colnames(dat)
gene = rownames(dat)
metric = cal_metric_Seurat(dat)
me[[1]] = metric

d.umap <- umap(t(dat))
df <- data.frame(x = d.umap$layout[,1],
                 y = d.umap$layout[,2],
                 Cell = cell)
ggplot(df, aes(x, y, colour = Cell)) +
  geom_point() +
  # ggtitle(methods_name) +
  labs(x="UMAP_1", y="UMAP_2") +
  theme_bw() +
  theme(
    axis.text.x=element_blank(),
    axis.text.y=element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(hjust = 0.5))

# scWMC
dat = readMat("D:/MyWorkWorld/Imputation/scWMC/scalability/Results/scWMC.mat")[[1]]
colnames(dat) = cell
rownames(dat) = gene
metric = cal_metric_Seurat(dat)
me[[7]] = metric

# ALRA
dat = readRDS(paste0(dirc, "ALRA.rds"))
colnames(dat) = cell
metric = cal_metric_Seurat(dat)
me[[2]] = metric

# DCA
dat = read.table(paste0(dirc, "DCA.tsv"), sep = "\t", header = T, row.names = 1)
colnames(dat) = cell
metric = cal_metric_Seurat(dat)
me[[3]] = metric

# MAGIC
dat = read.table(paste0(dirc, "MAGIC.csv"), sep = ",", header = T, row.names = 1)
colnames(dat) = cell
metric = cal_metric_Seurat(dat)
me[[4]] = metric

# SAVER
dat = readRDS(paste0(dirc, "SAVER.rds"))
colnames(dat) = cell
metric = cal_metric_Seurat(dat)
me[[5]] = metric

# scVI
dat = read.table(paste0(dirc, "scVI.csv"), sep = ",", header = F)
colnames(dat) = cell
metric = cal_metric_Seurat(dat)
me[[6]] = metric


saveRDS(me, file = "D:/MyWorkWorld/Imputation/scWMC/scalability/Results/NMI.rds")
