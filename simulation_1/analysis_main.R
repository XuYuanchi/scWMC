# This is the main R file to generate the results related our manuscript
# figure 2 and the supplementary figures related to figure 2.
# Please contact Tao Peng: pengt@email.chop.edu if you have any questions 
# about the scripts or data

# load the libraries 

library(dplyr)
library(vsn)
library(Seurat)
library(scater)
library(edgeR)
library(gridExtra)
library(R.matlab)
library(cowplot)
library(biomaRt)
library(data.table)
library(lattice)
library(scImpute)
library(SCRABBLE)
library(VennDiagram)
library(Rtsne)
library(DT)
library(ggpubr)
library(ggsignif)
library(scatterplot3d)
library(ggplot2)
library(reshape2)
library(ggfortify)
library(refGenome)
library(pheatmap)
library(RColorBrewer)
library(dendsort)
library(entropy)
library(DrImpute)
library(splatter)
library(RColorBrewer)
library(mcriPalettes)
library(plotly)
library(factoextra)
library(cluster)
library(NbClust)
library(fpc)
library(class)
library(VIPER)
library(SC3)
# load
library(Matrix)

source("analysis_library.R")


# the following script is to generate the simulation data. Here we use 
# HPC to generate the simulation which could reduce the running time
for(dropout_index in c(1:1) ){
  
  for(seed_value in c(1:1)){
    
    generate_save_data(dropout_index, seed_value)
    
  }
}


