library(scTSSR)

dat_name = c('Deng', 'Petropoulos')
for(i in c(2:2)){
  data = readRDS(paste0("D:/MyWorkWorld/Imputation/scWMC/ICT/data/", dat_name[i], "_filter.rds"))
  
  scTSSR_result <- scTSSR(data, percent=0.05, learning_rate=0.0001, epochs=100)
  saveRDS(scTSSR_result, file=paste0("D:/MyWorkWorld/Imputation/scWMC/ICT/Results/", dat_name[i], "/scTSSR.rds"))
}