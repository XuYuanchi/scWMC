library(DrImpute)

data_name = c("sc_CELseq2", "sc_10x", "sc_Droseq", "sc_10x_5cl", "sc_Celseq2_5cl_p1", "sc_Celseq2_5cl_p2", "sc_Celseq2_5cl_p3")
for (j in c(6:10)) {
  
for(i in c(3:6)){
    X = readRDS(paste0("D:/MyWorkWorld/Imputation/scWMC/clustering/data/", data_name[i], ".rds"))

    X <- log(X + 1)
    DrImpute_samp <- DrImpute(X)
    saveRDS(DrImpute_samp, file=paste0("D:/MyWorkWorld/Imputation/scWMC/Clustering/DrImpute/", j, "/", data_name[i], ".rds"))

}
}