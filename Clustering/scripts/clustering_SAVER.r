library(SAVER)

data_name = c("sc_CELseq2", "sc_10x", "sc_Droseq", "sc_10x_5cl", "sc_Celseq2_5cl_p1", "sc_Celseq2_5cl_p2", "sc_Celseq2_5cl_p3")
for (j in c(6:10)) {
    for (i in c(3:6)) {
        X = readRDS(paste0("/home/suyanchi/project/scWMC/Clustering/data/", data_name[i], ".rds"))

        saver_result <- saver(X, estimates.only = TRUE, ncores=30)
        saveRDS(saver_result, file = paste0("/home/suyanchi/project/scWMC/Clustering/SAVER/", j, "/", data_name[i], ".rds"))
    }
}
