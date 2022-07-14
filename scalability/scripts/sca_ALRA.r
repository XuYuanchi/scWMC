source('/home/grads/fuzhowang2/suyanchi/ALRA/alra.R')

data = as.matrix(read.table(file="/home/grads/fuzhowang2/suyanchi/scWMC/scalability/data/100k.csv", sep = ",", header = T, row.names=1))

A_norm <- normalize_data(t(data))
result.completed = alra(A_norm)[[3]]
ALRA_result = t(result.completed)
saveRDS(ALRA_result , "/home/grads/fuzhowang2/suyanchi/scWMC/scalability/Results/ALRA.rds")


