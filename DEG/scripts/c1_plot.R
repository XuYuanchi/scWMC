setwd("D:/MyWorkWorld/Imputation/scWMC/DEG/")
library(ggplot2)
library(ggsci)
library(RColorBrewer)
library(umap)
library(patchwork)
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(DESeq2)

color = readRDS("D:/MyWorkWorld/Imputation/scWMC/color.rds")
## F.score, AUC, ACC
load("Results/DEG.performance.RData")
load("Results/DEG.DESeq2.results.rd")
load("Results/DEG.DESeq2.scWMC.rd")
#load("Results/DEG.DESeq2.our_3_result.rd")
load("Results/DEG.DESeq2.MAGIC.rd")
load("Results/DEG.DESeq2.scTSSR.rd")
load("Results/DEG.DESeq2.ALRA.rd")
# load("data/DEG.DESeq2.DCA.rd")
load("Results/DEG.DESeq2.EnImpute.rd")
load("Results/DEG.DESeq2.PBLR.rd")
load("Results/DEG.DESeq2.scVI.rd")
set.seed(1234)
DEG.performance <- as.data.frame(DEG.performance)


methods <- c("Raw", "scWMC", "ALRA", "DCA", "DrImpute", "EnImpute", "MAGIC", "PBLR", "SAVER", 
             "scImpute", "scTSSR", "scVI", "VIPER")
methods <- factor(methods, levels = c("Raw", "scWMC", "ALRA", "DCA", "DrImpute", "EnImpute", "MAGIC", "PBLR", "SAVER", 
                                      "scImpute", "scTSSR", "scVI", "VIPER"))
DEG.performance$methods <- methods
p1 <- ggplot(DEG.performance, aes(x=methods, y=AUC)) + 
  coord_cartesian(ylim=c(0.8,0.9))+
  geom_bar(aes(fill = methods), stat = "identity") +
  scale_fill_manual(values = color) +
  labs( x="Methods", y="AUC") +
  theme_bw() +
  theme(
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    #axis.ticks.x=element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank()
  ) 
#scale_x_continuous(breaks=data$id) +
#guides(fill=FALSE)

p2 <- ggplot(DEG.performance, aes(x=methods, y=ACC)) + 
  geom_bar(aes(fill = methods), stat = "identity") +
  scale_fill_manual(values = color) +
  coord_cartesian(ylim=c(0.6,0.85))+
  labs(x="Methods",y="ACC") +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        #axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
#scale_x_continuous(breaks=data$id,labels=data$Methods) +
#guides(fill=FALSE)

p3 <- ggplot(DEG.performance, aes(x=methods, y=F.score)) + 
  geom_bar(aes(fill = methods), stat = "identity") +
  scale_fill_manual(values = color) +
  coord_cartesian(ylim=c(0,0.7))+
  labs(x="Methods",y="F.score") +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        #axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
#scale_x_continuous(breaks=data$id,labels=data$Methods) +
#guides(fill=FALSE)
#####################################################################

p = (p1 | p2 | p3) + plot_layout(guides = "collect") & theme(legend.position = "bottom", legend.key.size = unit(20, "pt")) 

saveRDS(p, file = "Results/metrix.rds")
#ggsave(filename = "plots/DEG_2.pdf", plot = p, width = 12, height = 4)


## up-down 1000 genes
plot_correlation <- function(object1, object2, xlab, ylab){
  ovlp <- intersect(rownames(object1[object1$padj <= 0.05,]),rownames(object2[object2$padj <= 0.05,]))
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

p4 <- list()
p4[[1]] <- plot_correlation(blk.result[blk.top1000.genes,], raw.result[rownames(raw.result) %in% blk.top1000.genes,], "Bulk", "Raw")
p4[[2]] <- plot_correlation(blk.result[blk.top1000.genes,], our.result[blk.top1000.genes,], "Bulk", "scWMC")
p4[[3]] <- plot_correlation(blk.result[blk.top1000.genes,], ALRA.result[intersect(ALRA.result@rownames, blk.top1000.genes),], "Bulk", "ALRA")
p4[[4]] <- plot_correlation(blk.result[blk.top1000.genes,], DCA.result[blk.top1000.genes,], "Bulk", "DCA")
p4[[5]] <- plot_correlation(blk.result[blk.top1000.genes,], DrImpute.result[blk.top1000.genes,], "Bulk", "DrImpute")
p4[[6]] <- plot_correlation(blk.result[blk.top1000.genes,], EnImpute.result[intersect(EnImpute.result@rownames, blk.top1000.genes),], "Bulk", "EnImpute")
p4[[7]] <- plot_correlation(blk.result[blk.top1000.genes,], MAGIC.result[rownames(MAGIC.result)%in%blk.top1000.genes,], "Bulk", "MAGIC")
p4[[8]] <- plot_correlation(blk.result[blk.top1000.genes,], PBLR.result[intersect(PBLR.result@rownames,blk.top1000.genes), ], "Bulk", "PBLR")
p4[[9]] <- plot_correlation(blk.result[blk.top1000.genes,], SAVER.result[rownames(SAVER.result) %in%blk.top1000.genes,], "Bulk", "SAVER")
p4[[10]] <- plot_correlation(blk.result[blk.top1000.genes,], scimpute.result[blk.top1000.genes,], "Bulk", "scImpute")
p4[[11]] <- plot_correlation(blk.result[blk.top1000.genes,], scTSSR.result[rownames(scTSSR.result)%in%blk.top1000.genes,], "Bulk", "scTSSR")                
p4[[12]]<- plot_correlation(blk.result[blk.top1000.genes,], scVI.result[blk.top1000.genes,], "Bulk", "scVI")
p4[[13]] <- plot_correlation(blk.result[blk.top1000.genes,], VIPER.result[blk.top1000.genes,], "Bulk", "VIPER")

layout1 <-'
ABCDEFGH
IJKLMNOP
'
p1 = p1 + theme(legend.position = "none")
p2 = p2 + theme(legend.position = "none")
p3 = p3 + theme(legend.position = "none")

p = p4[[1]] +  p4[[2]] + p4[[3]] + p4[[4]] + p4[[5]] + p4[[6]] + p4[[7]] + p4[[8]] + p4[[9]] + p4[[10]] + p4[[11]] + p4[[12]] + p4[[13]] + p1 +  p2 + p3 + 
  plot_layout(ncol = 8, guides = 'collect')  + plot_annotation(tag_levels = 'A')  & theme(legend.position = 'bottom')

saveRDS(p, file = "Results/DEG.rds")
saveRDS(p4, file = "Results/DEG_pearson.rds")
ggsave(filename = 'Results/DEG.pdf', plot = p, width = 16, height = 5)
