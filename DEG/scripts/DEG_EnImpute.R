# library(EnImpute)
source("/home/suyanchi/project/scWMC/utils/EnImpute.R")
source("/home/suyanchi/project/scWMC/utils/alra.R")
source("/home/suyanchi/project/scWMC/utils/magic.R")
source("/home/suyanchi/project/scWMC/utils/preprocessing.R")
source("/home/suyanchi/project/scWMC/utils//utils.R")

dat = read.table("/home/suyanchi/project/scWMC/DEG/data/DEG.raw.csv", sep = ",", header = F)

EnImpute_result = EnImpute(dat, ALRA = FALSE, DCA=FALSE, DrImpute=FALSE, MAGIC=FALSE, SAVER = FALSE, 
scImpute = TRUE, scRMD=FALSE, Seurat=FALSE, SAVER.ncores=20, scImpute.Kcluster=2)[2]

saveRDS(EnImpute_result, file="/home/suyanchi/project/scWMC/DEG/Results/EnImpute.rds")