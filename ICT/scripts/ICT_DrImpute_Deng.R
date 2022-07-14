library(DrImpute)

# read data

dat_name = c('Deng', 'Petropoulos')

for(i in c(2:2)){
  X = readRDS(paste0('D:/MyWorkWorld/Imputation/scWMC/ICT/data/', dat_name[i], '_filter.rds'))
  X <- preprocessSC(X, min.expressed.gene = 0)
  X <- log(X + 1)
  DrImpute_samp <- DrImpute(X)
  saveRDS(DrImpute_samp, file=paste0("D:/MyWorkWorld/Imputation/scWMC/ICT/Results/", dat_name[i], "/DrImpute.rds"))
  
}
