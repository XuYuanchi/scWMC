library(Seurat)
library(SingleCellExperiment)
library(SC3)
#library(mclust)
library(R.matlab)
library(aricode)
library(stats)
library(mclustcomp)
# define function
cal_metric <- function(count){
  
  sce = SingleCellExperiment(
    assays = list(
      counts = as.matrix(count),
      logcounts = log2(as.matrix(count) + 1)
    )
  )
  rowData(sce)$feature_symbol <- rownames(sce)
  k = length(unique(colnames(count)))
  kn=1000
  if(dim(count)[2]>2000){
    kn=50
  }
  sce = sc3(sce, ks = k, gene_filter=FALSE, kmeans_nstart=kn)
  km.label.k = sce@colData@listData[[1]]
  
  label <- colnames(count)
  
  #ARI <- adjustedRandIndex(as.factor(label), x.seurat$seurat_clusters)
  ari <- ARI(as.factor(label), km.label.k)
  nmi <- NMI(as.factor(label), km.label.k)
  jaccard <- mclustcomp(label, as.vector(km.label.k), types = "jaccard")$scores
  metric <- c(ari, nmi, jaccard)
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
  jaccard <- mclustcomp(as.vector(label), as.vector(x.seurat$seurat_clusters), types = "jaccard")$scores
  metric <- c(ari, nmi, jaccard)
  return (metric)
}
# data path
methods_name = c("raw", "scWMC", "ALRA", "DCA", "DrImpute", "EnImpute", "MAGIC", "PBLR", "SAVER", 
                 "scImpute", "scTSSR", "scVI", "VIPER")

dataset_name = c("sc_CELseq2", "sc_10x", "sc_Droseq", "sc_10x_5cl", "sc_Celseq2_5cl_p1", "sc_Celseq2_5cl_p2", "sc_Celseq2_5cl_p3")
cluster_numer = c(3, 3, 3, 5, 5, 5, 5)
dir = "D:/MyWorkWorld/Imputation/scWMC/Clustering/"
# ARI and NMI
ARI = matrix(nrow = 13, ncol = 7)
NMI = matrix(nrow = 13, ncol = 7)
jaccard = matrix(nrow = 13, ncol = 7)
rownames(ARI) = methods_name
colnames(ARI) = dataset_name
rownames(NMI) = methods_name
colnames(NMI) = dataset_name
rownames(jaccard) = methods_name
colnames(jaccard) = dataset_name
for (j in c(1:5)){
  for(i in c(3:6)) {
    raw_file = paste0(dir, "data/", dataset_name[i], ".rds")
    ALRA_file = paste0(dir, "ALRA/", j, "/", dataset_name[i], ".rds")
    DCA_file = paste0(dir, "DCA/", j, "/", dataset_name[i], ".tsv")
    DrImpute_file = paste0(dir, "DrImpute/", j, "/", dataset_name[i], ".rds")
    EnImpute_file = paste0(dir, "EnImpute/", j, "/", dataset_name[i], ".rds")
    MAGIC_file = paste0(dir, "MAGIC/", j, "/", dataset_name[i], ".csv")
    PBLR_file = paste0(dir, "PBLR/", j, "/", dataset_name[i], ".mat")
    SAVER_file = paste0(dir, "SAVER/", j, "/", dataset_name[i], ".rds")
    scImpute_file = paste0(dir, "scImpute/", j, "/", dataset_name[i], ".rds")
    scTSSR_file = paste0(dir, "scTSSR/", j, "/", dataset_name[i], ".rds")
    scVI_file = paste0(dir, "scVI/", j, "/", dataset_name[i], ".csv")
    VIPER_file = paste0(dir, "VIPER/", dataset_name[i], ".rds")
    scWMC_file = paste0(dir, "scWMC/", j, "/", dataset_name[i], "_m.mat")
    
    
    # raw data
    data = readRDS(raw_file)
    label = colnames(data)
    print(length(unique(label)))
    gene = rownames(data)
    # metric = cal_metric(data)
    # ARI[1, i] = metric[1]
    # NMI[1, i] = metric[2]
    # jaccard[1, i] = metric[3]
    # scWMC
    data = readMat(scWMC_file)[[1]]
    rownames(data) = gene
    colnames(data) = label
    metric = cal_metric(data)
    ARI[2, i] = metric[1]
    NMI[2, i] = metric[2]
    jaccard[2, i] = metric[3]
    # # ALRA
    # data = readRDS(ALRA_file)
    # colnames(data) = label
    # metric = cal_metric(data)
    # ARI[3, i] = metric[1]
    # NMI[3, i] = metric[2]
    # jaccard[3, i] = metric[3]
    # # DCA
    # data = read.table(DCA_file, header = T, sep = "\t")
    # rownames(data) = data[, 1]
    # data = data[, -1]
    # colnames(data) = label
    # metric = cal_metric(data)
    # ARI[4, i] = metric[1]
    # NMI[4, i] = metric[2]
    # jaccard[4, i] = metric[3]
    # # DrImpute
    # data = readRDS(DrImpute_file)
    # colnames(data) = label
    # metric = cal_metric(data)
    # ARI[5, i] = metric[1]
    # NMI[5, i] = metric[2]
    # jaccard[5, i] = metric[3]
    # # EnImpute
    # data = readRDS(EnImpute_file)[[1]]
    # colnames(data) = label
    # metric = cal_metric(data)
    # ARI[6, i] = metric[1]
    # NMI[6, i] = metric[2]
    # jaccard[6, i] = metric[3]
    # 
    # # MAGIC
    # data = read.table(MAGIC_file, header = TRUE, sep = ",")
    # data = data.matrix(data)[, -1]
    # rownames(data) = gene
    # colnames(data) = label
    # metric = cal_metric(data)
    # ARI[7, i] = metric[1]
    # NMI[7, i] = metric[2]
    # jaccard[7, i] = metric[3]
    # # PBLR
    # data = readMat(PBLR_file)[[1]]
    # rownames(data) = gene
    # colnames(data) = label
    # metric = cal_metric(data)
    # ARI[8, i] = metric[1]
    # NMI[8, i] = metric[2]
    # jaccard[8, i] = metric[3]
    # # SAVER
    # data = readRDS(SAVER_file)
    # colnames(data) = label
    # metric = cal_metric(data)
    # ARI[9, i] = metric[1]
    # NMI[9, i] = metric[2]
    # jaccard[9, i] = metric[3]
    # # scImpute
    # data = readRDS(scImpute_file)
    # colnames(data) = label
    # metric = cal_metric(data)
    # ARI[10, i] = metric[1]
    # NMI[10, i] = metric[2]
    # jaccard[10, i] = metric[3]
    # # scTSSR
    # data = readRDS(scTSSR_file)[[1]]
    # colnames(data) = label
    # metric = cal_metric(data)
    # ARI[11, i] = metric[1]
    # NMI[11, i] = metric[2]
    # jaccard[11, i] = metric[3]
    # # scVI
    # data = data.matrix(read.table(scVI_file, header = F, ","))
    # rownames(data) = gene
    # colnames(data) = label
    # metric = cal_metric(data)
    # ARI[12, i] = metric[1]
    # NMI[12, i] = metric[2]
    # jaccard[12, i] = metric[3]
    # # VIPER
    # data = readRDS(VIPER_file)[[1]]
    # colnames(data) = label
    # metric = cal_metric(data)
    # ARI[13, i] = metric[1]
    # NMI[13, i] = metric[2]
    # jaccard[13, i] = metric[3]
  }
  # save data
  saveRDS(ARI, file = paste0(dir, "Results/SC3/ARI_2", j, ".rds"))
  saveRDS(NMI, file = paste0(dir, "Results/SC3/NMI_2", j, ".rds"))
  saveRDS(jaccard, file = paste0(dir, "Results/SC3/jaccard_2", j, ".rds"))
}