#
source("D:/project/ALRA-master/alra.R")

# read data

dat_name = c('Deng', 'Petropoulos')

for(i in c(2:2)){
dat = readRDS(paste0('D:/MyWorkWorld/Imputation/scWMC/ICT/data/', dat_name[i], '_filter.rds'))

A_norm <- normalize_data(t(dat))
k_choice <- choose_k(A_norm)
A_norm_completed = alra(A_norm, k=k_choice$k)[[3]]
ALRA_result = t(A_norm_completed)
saveRDS(ALRA_result, file=paste0("D:/MyWorkWorld/Imputation/scWMC/ICT/Results/", dat_name[i], "/ALRA.rds"))
}