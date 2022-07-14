library(VIPER)

dat_name = c('Deng', 'Petropoulos')
for(i in c(2:2)){
  data = readRDS(paste0("D:/MyWorkWorld/Imputation/scWMC/ICT/data/", dat_name[i], "_filter.rds"))
  
  VIPER_samp <- VIPER(data, num = 1000, percentage.cutoff = 0.1, minbool = FALSE, 
                      alpha = 1, report = FALSE, outdir = NULL, prefix = NULL)
  saveRDS(VIPER_samp, file=paste0("D:/MyWorkWorld/Imputation/scWMC/ICT/Results/", dat_name[i], "/VIPER.rds"))
}