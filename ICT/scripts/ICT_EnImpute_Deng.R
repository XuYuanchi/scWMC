library(EnImpute)
#source("D:/project/EnImpute-master/pkg/R/EnImpute.R")
#source("D:/project/EnImpute-master/pkg/R/alra.R")

#source("D:/project/DrImpute-master/R/DrImpute.R")


dat_name = c('Deng', 'Petropoulos')

for(i in c(1:1)){
  dat = readRDS(paste0('/root/autodl-tmp/suyanchi/scWMC/ICT/data/', dat_name[i], '_raw.rds'))
  
  EnImpute_result = EnImpute(dat, ALRA = TRUE, DCA = FALSE, DrImpute = FALSE, MAGIC = FALSE, SAVER = FALSE, 
                             scImpute = FALSE, scRMD = FALSE, Seurat = FALSE, SAVER.ncores = 50, scImpute.ncores = 50)[2]
  
  saveRDS(EnImpute_result, file=paste0("/root/autodl-tmp/suyanchi/scWMC/ICT/Results/", dat_name[i], "/EnImpute.rds"))
  
}