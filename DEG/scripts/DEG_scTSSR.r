library(scTSSR)

X = read.table(file="/home/grads/fuzhowang2/suyanchi/scWMC/DEG/data/DEG.raw.tsv", sep="\t", header=T, row.names=1)

scTSSR_result <- scTSSR(X, percent=0, learning_rate=0.0001, epochs=10, estimates.only = TRUE, MAX_ITER = 1, ncores=12)

saveRDS(scTSSR_result, file="/home/grads/fuzhowang2/suyanchi/scWMC/DEG/Results/scTSSR.rds")
