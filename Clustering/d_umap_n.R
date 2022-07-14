library(umap)
library(ggplot2)
library(R.matlab)
library(patchwork)
set.seed(20)

plot_umap <- function(df, methods_name, Ja){
  p = ggplot(df, aes(x, y, colour = Cell)) +
    geom_point() +
    # ggtitle(methods_name) +
    labs(x="UMAP_1", y="UMAP_2") +
    ggtitle(methods_name, paste0("Jaccard Index = ", round(Ja, digits = 2)))
    theme_bw() +
    theme(
      axis.text.x=element_blank(),
      axis.text.y=element_blank(),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      plot.title = element_text(hjust = 0.5))
    
    return(p)
}


methods_name = c("Raw", "scWMC", "ALRA", "DCA", "DrImpute", "EnImpute", "MAGIC", "PBLR", "SAVER", 
                 "scImpute", "scTSSR", "scVI", "VIPER")

dataset_name = c("sc_CELseq2", "sc_10x", "sc_Droseq", "sc_10x_5cl", "sc_Celseq2_5cl_p1", "sc_Celseq2_5cl_p2", "sc_Celseq2_5cl_p3")
color = readRDS("D:/MyWorkWorld/Imputation/scWMC/color.rds")
dir = "D:/MyWorkWorld/Imputation/scWMC/Clustering/"
Ja = readRDS("D:/MyWorkWorld/Imputation/scWMC/Clustering/Re/Seurat/Ja_1.rds")


outdir = "plots/"
p = vector("list", 7)
for (j in c(5:5)){
  for(i in c(6:6)) {
    set.seed(20)
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
    scWMC_file = paste0(dir, "scWMC/", j,  "/", dataset_name[i], "_mm.mat")
    
    # raw data
    data = readRDS(raw_file)
    label = colnames(data)
    d.umap <- umap(t(data))
    df <- data.frame(x = d.umap$layout[,1],
                     y = d.umap$layout[,2],
                     Cell = label)
    
    p[[i]][[1]] = plot_umap(df, methods_name = methods_name[1], Ja[1, i])
    # scWMC
    data = readMat(scWMC_file)[[1]]
    colnames(data) = label
    d.umap <- umap(t(data))
    df <- data.frame(x = d.umap$layout[,1],
                     y = d.umap$layout[,2],
                     Cell = label)
    p[[i]][[2]] = plot_umap(df, methods_name = methods_name[2], Ja[2, i])
    
    # ALRA
    data = readRDS(ALRA_file)
    d.umap <- umap(t(data))
    df <- data.frame(x = d.umap$layout[,1],
                     y = d.umap$layout[,2],
                     Cell = label)
    
    p[[i]][[3]] = plot_umap(df, methods_name = methods_name[3], Ja[3, i])
    
    # DCA
    data = read.table(DCA_file, header = T, sep = "\t")
    rownames(data) = data[, 1]
    data = data[, -1]
    d.umap <- umap(t(data))
    df <- data.frame(x = d.umap$layout[,1],
                     y = d.umap$layout[,2],
                     Cell = label)
    
    p[[i]][[4]] = plot_umap(df, methods_name = methods_name[4], Ja[4, i])
    # DrImpute
    data = readRDS(DrImpute_file)
    colnames(data) = label
    d.umap <- umap(t(data))
    df <- data.frame(x = d.umap$layout[,1],
                     y = d.umap$layout[,2],
                     Cell = label)
    
    p[[i]][[5]] = plot_umap(df, methods_name = methods_name[5], Ja[5, i])
    # EnImpute
    data = readRDS(EnImpute_file)[[1]]
    d.umap <- umap(t(data))
    df <- data.frame(x = d.umap$layout[,1],
                     y = d.umap$layout[,2],
                     Cell = label)
    
    p[[i]][[6]] = plot_umap(df, methods_name = methods_name[6], Ja[6, i])
    # MAGIC
    data = read.table(MAGIC_file, header = TRUE, sep = ",", row.names = 1)
    d.umap <- umap(data)
    df <- data.frame(x = d.umap$layout[,1],
                     y = d.umap$layout[,2],
                     Cell = label)
    
    p[[i]][[7]] = plot_umap(df, methods_name = methods_name[7], Ja[7, i])
    # PBLR
    data = readMat(PBLR_file)[[1]]
    d.umap <- umap(t(data))
    df <- data.frame(x = d.umap$layout[,1],
                     y = d.umap$layout[,2],
                     Cell = label)
    
    p[[i]][[8]] = plot_umap(df, methods_name = methods_name[8], Ja[8, i])
    # SAVER
    data = readRDS(SAVER_file)
    d.umap <- umap(t(data))
    df <- data.frame(x = d.umap$layout[,1],
                     y = d.umap$layout[,2],
                     Cell = label)
    
    p[[i]][[9]] = plot_umap(df, methods_name = methods_name[9], Ja[9, i])
    # scImpute
    data = readRDS(scImpute_file)
    d.umap <- umap(t(data))
    df <- data.frame(x = d.umap$layout[,1],
                     y = d.umap$layout[,2],
                     Cell = label)
    
    p[[i]][[10]] = plot_umap(df, methods_name = methods_name[10], Ja[10, i])
    # scTSSR
    data = readRDS(scTSSR_file)[[1]]
    d.umap <- umap(t(data))
    df <- data.frame(x = d.umap$layout[,1],
                     y = d.umap$layout[,2],
                     Cell = label)
    
    p[[i]][[11]] = plot_umap(df, methods_name = methods_name[11], Ja[11, i])
    # scVI
    data = data.matrix(read.table(scVI_file, header = F, ","))
    d.umap <- umap(data)
    df <- data.frame(x = d.umap$layout[,1],
                     y = d.umap$layout[,2],
                     Cell = label)
    
    p[[i]][[12]] = plot_umap(df, methods_name = methods_name[12], Ja[12, i])
    # VIPER
    data = readRDS(VIPER_file)[[1]]
    d.umap <- umap(t(data))
    df <- data.frame(x = d.umap$layout[,1],
                     y = d.umap$layout[,2],
                     Cell = label)
    
    p[[i]][[13]] = plot_umap(df, methods_name = methods_name[13], Ja[13, i])
    
  }
}


layout <- "
ABCDEFG
HIJKLM#
"
p1 = p[[i]][[1]]
for (m in c(2:13)) {
  p1 = p1 + p[[i]][[m]]
}

p1 = p1 + plot_layout(guides = 'collect', design = layout) & theme(legend.position = 'bottom')
saveRDS(p1, file = paste0("D:/MyWorkWorld/Imputation/scWMC/Clustering/Re/Seurat/umap_seurat_", i, "_", j, ".rds"))
p1 = readRDS(paste0("D:/MyWorkWorld/Imputation/scWMC/Clustering/Re/Seurat/umap_seurat_", i, "_", j, ".rds"))
ggsave(filename = paste0("D:/MyWorkWorld/Imputation/scWMC/Clustering/Re/Seurat/umap_seurat_", i, "_", j, ".pdf"), plot = p1, width = 14, height = 5)


p2 = readRDS("D:/MyWorkWorld/Imputation/scWMC/Clustering/Re/NMI_Seurat.rds")

design = "
ABCD
EEEE"

p3 = p2[[1]][[3]] + p2[[1]][[4]] + p2[[1]][[5]] + p2[[1]][[6]] + p1  +
  plot_layout(guides = "collect", design=design)  & theme(legend.position = "bottom")
ggsave(filename = paste0("D:/MyWorkWorld/Imputation/scWMC/Clustering/Re/Seurat/seurat_", i, "_", j, ".pdf"), plot = p3, width = 16, height = 9)
