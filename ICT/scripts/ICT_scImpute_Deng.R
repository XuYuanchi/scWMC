library(scImpute)

dat_name = c('Deng', 'Petropoulos')
for(i in c(1:1)){
  filename = paste0('D:/MyWorkWorld/Imputation/scWMC/ICT/data/', dat_name[i], '.csv')
  
  scimpute(count_path = filename, infile = "csv", outfile = "rds", out_dir = "./", 
           labeled = FALSE, drop_thre = 0.5, Kcluster = 10, ncores = 1)
  
}