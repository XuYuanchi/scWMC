library(scImpute)

data_name = c("sc_CELseq2", "sc_10x", "sc_Droseq", "sc_10x_5cl", "sc_Celseq2_5cl_p1", "sc_Celseq2_5cl_p2", "sc_Celseq2_5cl_p3")
for (j in c(1:10)) {
    for (i in c(3:6)) {
        filename = paste0("/home/suyanchi/project/scWMC/Clustering/data/", data_name[i], ".csv")
        scimpute(count_path = filename, infile = "csv", outfile = "rds", out_dir = paste0("/home/suyanchi/project/scWMC/Clustering/scImpute/", j, "/", data_name[i], "/"), 
labeled = FALSE, drop_thre = 0.5, Kcluster = 5, ncores = 5)
        dat = readRDS(paste0("/home/suyanchi/project/scWMC/Clustering/scImpute/", j, "/", data_name[i], "/scimpute_count.rds"))
        saveRDS(dat, file=paste0("/home/suyanchi/project/scWMC/Clustering/scImpute/", j, "/", data_name[i], ".rds"))
    }
}