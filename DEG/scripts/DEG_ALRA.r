#
source("/home/grads/fuzhowang2/suyanchi/ALRA/alra.R")

# read data

dat = read.table("/home/grads/fuzhowang2/suyanchi/scWMC/DEG/data/DEG.raw.csv", sep = ",", header = F)

A_norm <- normalize_data(t(dat))
k_choice <- choose_k(A_norm)
A_norm_completed = alra(A_norm, k=k_choice$k)[[3]]
ALRA_result = t(A_norm_completed)
saveRDS(ALRA_result, file="/home/grads/fuzhowang2/suyanchi/scWMC/DEG/Results/ALRA.rds")