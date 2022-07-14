library(EnImpute)
#source("/home/suyanchi/project/scWMC/utils/EnImpute.R")
#source("/home/suyanchi/project/scWMC/utils/alra.R")
#source("/home/suyanchi/project/scWMC/utils/magic.R")
#source("/home/suyanchi/project/scWMC/utils/preprocessing.R")
#source("/home/suyanchi/project/scWMC/utils//utils.R")

data_name = c("sc_CELseq2", "sc_10x", "sc_Droseq", "sc_10x_5cl", "sc_Celseq2_5cl_p1", "sc_Celseq2_5cl_p2", "sc_Celseq2_5cl_p3")
for(j in c(1:10)){
for(i in c(4)){
  data = readRDS(paste0("/home/suyanchi/project/scWMC/Clustering/data/", data_name[i], ".rds"))
  
  EnImpute_result = EnImpute(data, Seurat=FALSE)[2]
  saveRDS(EnImpute_result, file=paste0("/home/suyanchi/project/scWMC/Clustering/EnImpute/", j, "/", data_name[i], ".rds"))
}}