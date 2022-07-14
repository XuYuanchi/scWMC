library(VIPER)


data_name = c("sc_CELseq2", "sc_10x", "sc_Droseq", "sc_10x_5cl", "sc_Celseq2_5cl_p1", "sc_Celseq2_5cl_p2", "sc_Celseq2_5cl_p3")
for(j in c(7:10)){
    for(i in c(4:4)){
    data = readRDS(paste0("D:/MyWorkWorld/Imputation/scWMC/clustering/data/", data_name[i], ".rds"))

    VIPER_samp <- VIPER(data, num = 5000, percentage.cutoff = 0.1, minbool = FALSE, alpha = 1, report = FALSE, outdir = NULL, prefix = NULL)
    saveRDS(VIPER_samp, file=paste0("D:/MyWorkWorld/Imputation/scWMC/clustering/VIPER/", j, "/", data_name[i], ".rds"))

}
}