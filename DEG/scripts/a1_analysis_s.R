library(DESeq2)
library(R.matlab)

plot_correlation <- function(object1, object2, xlab, ylab){
  ovlp <- intersect(rownames(object1[object1$padj <= 0.05,]),rownames(object2[!is.na(object2$padj) &object2$padj <= 0.05,]))
  x = object1[ovlp,]$log2FoldChange
  y = object2[ovlp,]$log2FoldChange
  cor <- cor.test(x,y)
  r <- round(cor$estimate,2) 
  p <- signif(cor$p.value,3)
  up <- sum(x>0)
  down <- sum(x<0)
  data <- data.frame(x=x, y=y)
  p <- data %>% mutate(Color = ifelse(x*y>0, ifelse(x+y > 0, "#E64B35FF", "#3C5488FF"), "#F39B7FFF")) %>%
    ggplot(aes(x=x, y=y, color = Color)) +
    geom_point() +
    stat_smooth(method = 'lm', colour = 'orange') +
    annotate("text", x = -Inf, y = Inf, label = paste0("r = ", r), hjust = -0.2, 
             vjust = 2, size = 3) +
    annotate("text", x = -Inf, y = Inf, label = paste0("up = ", up), hjust = -0.2, 
             vjust = 4, size = 3) +
    annotate("text", x = -Inf, y = Inf, label = paste0("down = ", down), hjust = -0.2, 
             vjust = 6, size = 3) +
    scale_color_identity() +
    xlab(NULL) +
    theme_bw() +
    labs(x = xlab, y = ylab) +
    theme(
      axis.title.y = element_text(size = 12),
      axis.text = element_text(size = 12),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_blank()
    )
  return(p)
}

load("Results/DEG.DESeq2.results.rd")



setwd('D:/MyWorkWorld/Imputation/scWMC/DEG/')
dir = "D:/MyWorkWorld/Imputation/scWMC/DEG/"
label.file = paste0(dir, "data/","DEG.label.txt")
# raw 
raw.cts = read.table(paste0(dir, "data/DEG.raw.tsv"), sep = "\t", stringsAsFactors = F)
raw.cts <- round(raw.cts)
raw.label <- read.table(label.file)
row.names(raw.label) <- colnames(raw.cts)
colnames(raw.label) <- "Cell"

#scWMC
our.cts <- readMat('Results/scWMC.mat')[[1]]
rownames(our.cts) <- rownames(raw.cts)
colnames(our.cts) <- colnames(raw.cts)
our.cts <- round(data.matrix(our.cts))
our.dds <- DESeqDataSetFromMatrix(countData = our.cts,
                                  colData = raw.label,
                                  design = ~ Cell)

our.deg <- DESeq(our.dds)
our.result <- results(our.deg, contrast = c('Cell', 'H1', 'DEC'))
our.result <- our.result[complete.cases(our.result),]
our.deg.list <- our.result[our.result$padj <= 0.05 & abs(our.result$log2FoldChange) >1.5 & our.result$baseMean >= 10,]
our.deg.up.num <- sum(our.deg.list$log2FoldChange>0)
our.deg.down.num <- sum(our.deg.list$log2FoldChange<0)

our.deg.genes <- rownames(our.deg.list)
save(our.deg.genes, our.result, our.deg, file = "Results/DEG.DESeq2.scWMC.rd")

plot_correlation(blk.result[blk.top1000.genes,], our.result[blk.top1000.genes,], "Bulk", "scWMC")



## EnImpute
EnImpute.cts <- readRDS(paste0(dir, 'Results/EnImpute.rds'))[[1]]
rownames(EnImpute.cts) <- rownames(raw.cts)
colnames(EnImpute.cts) <- colnames(raw.cts)
EnImpute.cts <- round(data.matrix(EnImpute.cts))
EnImpute.cts <- DESeqDataSetFromMatrix(countData = EnImpute.cts,
                                     colData = raw.label,
                                     design = ~ Cell)

EnImpute.deg <- DESeq(EnImpute.cts)
EnImpute.result <- results(EnImpute.deg, contrast = c('Cell', 'H1', 'DEC'))
EnImpute.result <- EnImpute.result[complete.cases(EnImpute.result),]
EnImpute.deg.list <- EnImpute.result[EnImpute.result$padj <= 0.05 & abs(EnImpute.result$log2FoldChange) >1.5 & EnImpute.result$baseMean >= 10,]
EnImpute.deg.up.num <- sum(EnImpute.deg.list$log2FoldChange>0)
EnImpute.deg.down.num <- sum(EnImpute.deg.list$log2FoldChange<0)

EnImpute.deg.genes <- rownames(EnImpute.deg.list)
save(EnImpute.deg.genes, EnImpute.result, EnImpute.deg, file = paste0(dir, "Results/DEG.DESeq2.EnImpute.rd"))

## PBLR
PBLR.cts <- readMat(paste0(dir, 'Results/PBLR.mat'))[[1]]
rownames(PBLR.cts) <- rownames(raw.cts)
colnames(PBLR.cts) <- colnames(raw.cts)
PBLR.cts <- round(data.matrix(10^PBLR.cts-1))
PBLR.cts <- DESeqDataSetFromMatrix(countData = PBLR.cts,
                                       colData = raw.label,
                                       design = ~ Cell)

PBLR.deg <- DESeq(PBLR.cts)
PBLR.result <- results(PBLR.deg, contrast = c('Cell', 'H1', 'DEC'))
PBLR.result <- PBLR.result[complete.cases(PBLR.result),]
PBLR.deg.list <- PBLR.result[PBLR.result$padj <= 0.05 & abs(PBLR.result$log2FoldChange) >1.5 & PBLR.result$baseMean >= 10,]
PBLR.deg.up.num <- sum(PBLR.deg.list$log2FoldChange>0)
PBLR.deg.down.num <- sum(PBLR.deg.list$log2FoldChange<0)

PBLR.deg.genes <- rownames(PBLR.deg.list)
save(PBLR.deg.genes, PBLR.result, PBLR.deg, file = paste0(dir, "Results/DEG.DESeq2.PBLR.rd"))


## scVI
scVI.cts <- read.table(paste0(dir, 'Results/scVI.csv'), sep = ',', header = F)
rownames(scVI.cts) <- rownames(raw.cts)
colnames(scVI.cts) <- colnames(raw.cts)
scVI.cts <- round(data.matrix(10000*scVI.cts))
scVI.cts <- DESeqDataSetFromMatrix(countData = scVI.cts,
                                   colData = raw.label,
                                   design = ~ Cell)

scVI.deg <- DESeq(scVI.cts)
scVI.result <- results(scVI.deg, contrast = c('Cell', 'H1', 'DEC'))
scVI.result <- scVI.result[complete.cases(scVI.result),]
scVI.deg.list <- scVI.result[scVI.result$padj <= 0.05 & abs(scVI.result$log2FoldChange) >1.5 & scVI.result$baseMean >= 10,]
scVI.deg.up.num <- sum(scVI.deg.list$log2FoldChange>0)
scVI.deg.down.num <- sum(scVI.deg.list$log2FoldChange<0)

scVI.deg.genes <- rownames(scVI.deg.list)
save(scVI.deg.genes, scVI.result, scVI.deg, file = paste0(dir, "Results/DEG.DESeq2.scVI.rd"))

