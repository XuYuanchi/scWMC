## defined a function to compute the precision, recall, and F value of two cluster results.
setwd("D:/MyWorkWorld/Imputation/scWMC/DEG/")
library(ROCR)
precision.recall <- function(c1, c2, plot = F, label="ROC.curve"){
  #ref.label <- as.vector(outer(c1, c1, "=="))
  #pred.label <- as.vector(outer(c2, c2, "=="))
  #pred <- prediction(as.numeric(pred.label), as.numeric(ref.label) )
  pred <- prediction(as.numeric(c1), as.numeric(c2))
  
  # Accuracy            
  acc.tmp <- performance(pred, "acc");
  acc <- as.numeric(acc.tmp@y.values[[1]][2])
  
  # ROC curve
  ROC.perf <- performance(pred, "tpr", "fpr");
  # ROC area under the curve
  auc.tmp <- performance(pred,"auc");
  auc <- as.numeric(auc.tmp@y.values)
  
  if(plot){
    #pdf(paste(label,".pdf", sep = ""),2.5, 3.0)
    plot (ROC.perf);
    text(0.6, 0.05, paste("AUC=",round(auc, digits=4)))
    abline(a=0, b= 1)
    #dev.off()
  }
  
  ##precision
  prec.tmp = performance(pred, "ppv")
  prec <- as.numeric(prec.tmp@y.values[[1]][2])
  ## F1-score
  f.tmp = performance(pred, "f")
  f <- as.numeric(f.tmp@y.values[[1]][2])
  ## 
  return(list(F.score=f, AUC=auc, ACC=acc))
}

load("Results/DEG.DESeq2.results.rd")
# load("data/DEG.DESeq2.our_scimpute_result.rd")
load("Results/DEG.DESeq2.scWMC.rd")
load("Results/DEG.DESeq2.scTSSR.rd")
load("Results/DEG.DESeq2.MAGIC.rd")
load("Results/DEG.DESeq2.ALRA.rd")
load("Results/DEG.DESeq2.EnImpute.rd")
load("Results/DEG.DESeq2.PBLR.rd")
load("Results/DEG.DESeq2.scVI.rd")
# 
lfc=1.5

blk.deg.genes = rownames(blk.result[blk.result$padj <= 0.05 & abs(blk.result$log2FoldChange) >lfc & blk.result$baseMean >= 10, ])
our.deg.genes= rownames(our.result[our.result$padj <= 0.05 & abs(our.result$log2FoldChange) >lfc & our.result$baseMean >= 10, ])
ALRA.deg.genes= rownames(ALRA.result[ALRA.result$padj <= 0.05 & abs(ALRA.result$log2FoldChange) >lfc & ALRA.result$baseMean >= 10, ])
DCA.deg.genes= rownames(DCA.result[DCA.result$padj <= 0.05 & abs(DCA.result$log2FoldChange) >lfc & DCA.result$baseMean >= 10, ])
DrImpute.deg.genes= rownames(DrImpute.result[DrImpute.result$padj <= 0.05 & abs(DrImpute.result$log2FoldChange) >lfc & DrImpute.result$baseMean >= 10, ])
EnImpute.deg.genes= rownames(EnImpute.result[EnImpute.result$padj <= 0.05 & abs(EnImpute.result$log2FoldChange) >lfc & EnImpute.result$baseMean >= 10, ])
MAGIC.deg.genes= rownames(MAGIC.result[MAGIC.result$padj <= 0.05 & abs(MAGIC.result$log2FoldChange) >lfc & MAGIC.result$baseMean >= 10, ])
PBLR.deg.genes= rownames(PBLR.result[PBLR.result$padj <= 0.05 & abs(PBLR.result$log2FoldChange) >lfc & PBLR.result$baseMean >= 10, ])
SAVER.deg.genes= rownames(SAVER.result[SAVER.result$padj <= 0.05 & abs(SAVER.result$log2FoldChange) >lfc & SAVER.result$baseMean >= 10, ])
scimpute.deg.genes= rownames(scimpute.result[scimpute.result$padj <= 0.05 & abs(scimpute.result$log2FoldChange) >lfc & scimpute.result$baseMean >= 10, ])
scTSSR.deg.genes= rownames(scTSSR.result[scTSSR.result$padj <= 0.05 & abs(scTSSR.result$log2FoldChange) >lfc & scTSSR.result$baseMean >= 10, ])
scVI.deg.genes= rownames(scVI.result[scVI.result$padj <= 0.05 & abs(scVI.result$log2FoldChange) >lfc & scVI.result$baseMean >= 10, ])
VIPER.deg.genes= rownames(VIPER.result[VIPER.result$padj <= 0.05 & abs(VIPER.result$log2FoldChange) >lfc & VIPER.result$baseMean >= 10, ])

# p_v = 0.05
# 
# blk.deg.genes = rownames(blk.result[blk.result$padj <= p_v , ])
# our.deg.genes= rownames(our.result[our.result$padj <= p_v, ])
# ALRA.deg.genes= rownames(ALRA.result[ALRA.result$padj <= p_v , ])
# DCA.deg.genes= rownames(DCA.result[DCA.result$padj <= p_v , ])
# DrImpute.deg.genes= rownames(DrImpute.result[DrImpute.result$padj <= p_v, ])
# EnImpute.deg.genes= rownames(EnImpute.result[EnImpute.result$padj <= p_v, ])
# MAGIC.deg.genes= rownames(MAGIC.result[MAGIC.result$padj <= p_v , ])
# PBLR.deg.genes= rownames(PBLR.result[PBLR.result$padj <= p_v , ])
# SAVER.deg.genes= rownames(SAVER.result[SAVER.result$padj <= p_v, ])
# scimpute.deg.genes= rownames(scimpute.result[scimpute.result$padj <= p_v , ])
# scTSSR.deg.genes= rownames(scTSSR.result[scTSSR.result$padj <= p_v, ])
# scVI.deg.genes= rownames(scVI.result[scVI.result$padj <= p_v , ])
# VIPER.deg.genes= rownames(VIPER.result[VIPER.result$padj <= p_v, ])

all.degs = unique(c(blk.deg.genes , 
                    raw.deg.genes , 
                    DeepImpute.deg.genes, 
                    DrImpute.deg.genes,
                    scimpute.deg.genes, 
                    SAVER.deg.genes, 
                    DCA.deg.genes, 
                    VIPER.deg.genes,
                    scTSSR.deg.genes,
                    MAGIC.deg.genes,
                    ALRA.deg.genes,
                    EnImpute.deg.genes,
                    PBLR.deg.genes,
                    scVI.deg.genes,
                    our.deg.genes))

all.gene = rownames(blk.cts)

blk.deg.pred = all.gene %in% all.degs
raw.deg.pred = all.gene %in% raw.deg.genes
our.deg.pred = all.gene %in% our.deg.genes
ALRA.deg.pred = all.gene %in% ALRA.deg.genes
DCA.deg.pred = all.gene %in% DCA.deg.genes
DrImpute.deg.pred = all.gene %in% DrImpute.deg.genes
EnImpute.deg.pred = all.gene %in% EnImpute.deg.genes
MAGIC.deg.pred = all.gene %in% MAGIC.deg.genes
PBLR.deg.pred = all.gene %in% PBLR.deg.genes
SAVER.deg.pred = all.gene %in% SAVER.deg.genes
scimpute.deg.pred = all.gene %in% scimpute.deg.genes
scTSSR.deg.pred = all.gene %in% scTSSR.deg.genes
scVI.deg.pred = all.gene %in% scVI.deg.genes
VIPER.deg.pred = all.gene %in% VIPER.deg.genes


DEG.performance = rbind(
  unlist(precision.recall(blk.deg.pred, raw.deg.pred)),
  unlist(precision.recall(blk.deg.pred, our.deg.pred)),
  unlist(precision.recall(blk.deg.pred, ALRA.deg.pred)),
  unlist(precision.recall(blk.deg.pred, DCA.deg.pred)),
  unlist(precision.recall(blk.deg.pred, DrImpute.deg.pred)),
  unlist(precision.recall(blk.deg.pred, EnImpute.deg.pred)),
  unlist(precision.recall(blk.deg.pred, MAGIC.deg.pred)),
  unlist(precision.recall(blk.deg.pred, PBLR.deg.pred)),
  unlist(precision.recall(blk.deg.pred, SAVER.deg.pred)),
  unlist(precision.recall(blk.deg.pred, scimpute.deg.pred)),
  unlist(precision.recall(blk.deg.pred, scTSSR.deg.pred)),
  unlist(precision.recall(blk.deg.pred, scVI.deg.pred)),
  unlist(precision.recall(blk.deg.pred, VIPER.deg.pred))
)
rownames(DEG.performance) <- c("Raw", "scWMC", "ALRA", "DCA", "DrImpute", "EnImpute", "MAGIC", "PBLR", "SAVER", 
                               "scImpute", "scTSSR", "scVI", "VIPER")

save(DEG.performance, file="Results/DEG.performance.RData")
