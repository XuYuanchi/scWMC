library(SAVER)

dat_name = c('Deng', 'Petropoulos')
for(i in c(2:2)){
  data = readRDS(paste0('D:/MyWorkWorld/Imputation/scWMC/ICT/data/', dat_name[i], '_filter.rds'))
  
  saver_result <- saver(data, estimates.only = TRUE, ncores = 6)
  saveRDS(saver_result, file=paste0("D:/MyWorkWorld/Imputation/scWMC/ICT/Results/", dat_name[i], "/SAVER.rds"))
}