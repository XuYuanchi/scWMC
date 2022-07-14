source("D:/MyWorkWorld/Imputation/ALRA/alra.R")

data_name = c("sc_CELseq2", "sc_10x", "sc_Droseq", "sc_10x_5cl", "sc_Celseq2_5cl_p1", "sc_Celseq2_5cl_p2", "sc_Celseq2_5cl_p3")

for (j in c(6:10)) {
  

for (i in c(3:6)){
  data = readRDS(paste0("D:/MyWorkWorld/Imputation/scWMC/clustering/data/", data_name[i], ".rds"))
  A_norm <- normalize_data(t(data))
  k_choice <- choose_k(A_norm)
  A_norm_completed = alra(A_norm, k=k_choice$k)[[3]]
  ALRA_result = t(A_norm_completed)
  saveRDS(ALRA_result, file=paste0("D:/MyWorkWorld/Imputation/scWMC/clustering/ALRA/", j, "/", data_name[i], ".rds"))
}
}