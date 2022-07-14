library(scTSSR)

data_name = c("sc_CELseq2", "sc_10x", "sc_Droseq", "sc_10x_5cl", "sc_Celseq2_5cl_p1", "sc_Celseq2_5cl_p2", "sc_Celseq2_5cl_p3")
for (j in c(6:10)) {
    for (i in c(3:6)) {
        X = readRDS(paste0("/home/suyanchi/project/scWMC/Clustering/data/", data_name[i], ".rds"))

        scTSSR_result <- scTSSR(X, percent = 0.05, learning_rate = 0.0001, epochs = 100)
        saveRDS(scTSSR_result, file = paste0("/home/suyanchi/project/scWMC/Clustering/scTSSR/", j, "/", data_name[i], ".rds"))
    }
}