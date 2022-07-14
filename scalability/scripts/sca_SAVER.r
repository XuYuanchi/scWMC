library(SAVER)

data = as.matrix(read.table(file="/home/grads/fuzhowang2/suyanchi/scWMC/scalability/data/100k.csv", sep = ",", header = T, row.names=1))
saver_result <- saver(data, ncores = 20, estimates.only = TRUE)
saveRDS(saver_result , "/home/grads/fuzhowang2/suyanchi/scWMC/scalability/Results/SAVER.rds")