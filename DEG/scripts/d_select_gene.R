setwd("D:/MyWorkWorld/Imputation/scWMC/DEG/")

load("Results/DEG.DESeq2.results.rd")
# load("data/DEG.DESeq2.our_scimpute_result.rd")
load("Results/DEG.DESeq2.scWMC.rd")
load("Results/DEG.DESeq2.scTSSR.rd")
load("Results/DEG.DESeq2.MAGIC.rd")
load("Results/DEG.DESeq2.ALRA.rd")
load("Results/DEG.DESeq2.EnImpute.rd")
load("Results/DEG.DESeq2.PBLR.rd")
load("Results/DEG.DESeq2.scVI.rd")

# select DE
# 
padj = 0.01
lfc=1.5
basemean = 10

blk.deg.list = blk.result[blk.result$padj <= padj & abs(blk.result$log2FoldChange) >lfc & blk.result$baseMean >= 10, ]
raw.deg.list = raw.result[raw.result$padj <= padj & abs(raw.result$log2FoldChange) >lfc & raw.result$baseMean >= 10, ]
our.deg.list= our.result[our.result$padj <= padj & abs(our.result$log2FoldChange) >lfc & our.result$baseMean >= 10, ]
ALRA.deg.list= ALRA.result[ALRA.result$padj <= padj & abs(ALRA.result$log2FoldChange) >lfc & ALRA.result$baseMean >= 10, ]
DCA.deg.list= DCA.result[DCA.result$padj <= padj & abs(DCA.result$log2FoldChange) >lfc & DCA.result$baseMean >= 10, ]
DrImpute.deg.list= DrImpute.result[DrImpute.result$padj <= padj & abs(DrImpute.result$log2FoldChange) >lfc & DrImpute.result$baseMean >= 10, ]
EnImpute.deg.list= EnImpute.result[EnImpute.result$padj <= padj & abs(EnImpute.result$log2FoldChange) >lfc & EnImpute.result$baseMean >= 10, ]
MAGIC.deg.list= MAGIC.result[MAGIC.result$padj <= padj & abs(MAGIC.result$log2FoldChange) >lfc & MAGIC.result$baseMean >= 10, ]
PBLR.deg.list= PBLR.result[PBLR.result$padj <= padj & abs(PBLR.result$log2FoldChange) >lfc & PBLR.result$baseMean >= 10, ]
SAVER.deg.list= SAVER.result[SAVER.result$padj <= padj & abs(SAVER.result$log2FoldChange) >lfc & SAVER.result$baseMean >= 10, ]
scimpute.deg.list= scimpute.result[scimpute.result$padj <= padj & abs(scimpute.result$log2FoldChange) >lfc & scimpute.result$baseMean >= 10, ]
scTSSR.deg.list= scTSSR.result[scTSSR.result$padj <= padj & abs(scTSSR.result$log2FoldChange) >lfc & scTSSR.result$baseMean >= 10, ]
scVI.deg.list= scVI.result[scVI.result$padj <= padj & abs(scVI.result$log2FoldChange) >lfc & scVI.result$baseMean >= 10, ]
VIPER.deg.list= VIPER.result[VIPER.result$padj <= padj & abs(VIPER.result$log2FoldChange) >lfc & VIPER.result$baseMean >= 10, ]

blk_Ordered <- blk.deg.list[order(blk.deg.list$padj),]
raw_Ordered <- raw.deg.list[order(raw.deg.list),]
scWMC_Ordered <- our.deg.list[order(our.deg.list$padj),]
ALRA_Ordered <- ALRA.deg.list[order(ALRA.deg.list$padj),]
DCA_Ordered <- DCA.deg.list[order(DCA.deg.list$padj),]
DrImpute_Ordered <- DrImpute.deg.list[order(DrImpute.deg.list$padj),]
EnImpute_Ordered <- EnImpute.deg.list[order(EnImpute.deg.list$padj),]
MAGIC_Ordered <- MAGIC.deg.list[order(MAGIC.deg.list$padj),]
PBLR_Ordered <- PBLR.deg.list[order(PBLR.deg.list$padj),]
SAVER_Ordered <- SAVER.deg.list[order(SAVER.deg.list$padj),]
scimpute_Ordered <- scimpute.deg.list[order(scimpute.deg.list$padj),]
scTSSR_Ordered <- scTSSR.deg.list[order(scTSSR.deg.list$padj),]
scVI_Ordered <- scVI.deg.list[order(scVI.deg.list$padj),]
VIPER_Ordered <- VIPER.deg.list[order(VIPER.deg.list$padj),]

gene_name = list()
gene_name[[1]] = rownames(blk_Ordered)
gene_name[[2]] = rownames(raw_Ordered)
gene_name[[3]] = rownames(scWMC_Ordered)
gene_name[[4]] = rownames(ALRA_Ordered)
gene_name[[5]] = rownames(DCA_Ordered)
gene_name[[6]] = rownames(DrImpute_Ordered)
gene_name[[7]] = rownames(EnImpute_Ordered)
gene_name[[8]] = rownames(MAGIC_Ordered)
gene_name[[9]] = rownames(PBLR_Ordered)
gene_name[[10]] = rownames(SAVER_Ordered)
gene_name[[11]] = rownames(scimpute_Ordered)
gene_name[[12]] = rownames(scTSSR_Ordered)
gene_name[[13]] = rownames(scVI_Ordered)
gene_name[[14]] = rownames(VIPER_Ordered)

saveRDS(gene_name, file = "Results/gene_name.rds")

overlap_num = matrix(nrow = 10,  ncol = 13)

colnames(overlap_num) = c("Raw", "scWMC", "ALRA", "DCA", "DrImpute", "EnImpute", "MAGIC", "PBLR", "SAVER", 
                          "scImpute", "scTSSR", "scVI", "VIPER")
rownames(overlap_num) = c(1:10)


for (i in c(1:6)) {
  
  n = 50 * i
  
  for (j in c(2:14)) {
    
    overlap_num[i, j-1] = length(intersect(gene_name[[j]][1:n], gene_name[[1]][1:50]))
    
  }
  
}

saveRDS(overlap_num, file = "Results/blk50.rds")


