setwd("D:/MyWorkWorld/Imputation/scWMC/simulation_1/")

library(R.matlab)


source("a_calculate_library.r")

for(dropout_index in c(1:3) ){
  
  for(seed_value in c(1:10)){
    
    result <- run_error(dropout_index, seed_value)
    # result <- calculate_error_splatter(drop_index, seed_value)
    
    dir.create(file.path('error_data'), showWarnings = FALSE)
    
    saveRDS(result,
            file = paste0("error_data/error_",dropout_index,"_",seed_value,".rds")
    )
    
  }
  
}

# the following script is to calculate the errors. Here we use 
# HPC to impute the data using scrabble which could reduce the running time
for(dropout_index in c(1:3) ){
  
  for(seed_value in c(1:10)){
    
    cal_cell_distribution(dropout_index, seed_value)
    
  }
  
}
# the following script is to calculate the errors. Here we use 
# HPC to impute the data using scrabble which could reduce the running time
for(dropout_index in c(1:3) ){
  
  for(seed_value in c(1:10)){
    
    cal_gene_distribution(dropout_index, seed_value)
    
  }
  
}
# Gather the errors
# -----------------
error_list <- list()

error_cell_list <- list()

error_gene_list <- list()

# gather the error from the data from different dropout rates
for(i in c(1:3)){
  error_matrix <- c()
  
  error_cell_matrix <- c()
  
  error_gene_matrix <- c()
  
  for(j in c(1:10)){
    
    tmp <- readRDS(file = paste0("error_data/error_",i,"_",j,".rds"))
    
    error_matrix <- cbind(error_matrix,as.matrix(tmp$error))
    
    error_cell_matrix <- cbind(error_cell_matrix,as.matrix(tmp$error_cell))
    
    error_gene_matrix <- cbind(error_gene_matrix,as.matrix(tmp$error_gene))
    
  }
  
  error_list[[i]] <- error_matrix
  
  error_cell_list[[i]] <- error_cell_matrix
  
  error_gene_list[[i]] <- error_gene_matrix
  
}

# save the errors
saveRDS(error_list,file = "data/error_all.rds")

saveRDS(error_cell_list,file = "data/error_all_cell.rds")

saveRDS(error_gene_list,file = "data/error_all_gene.rds")
# ------

gene_error_list = list()
gene_error_cell_list = list()
gene_error_gene_list = list()
gene_error = matrix(nrow = 5, ncol = 100)
gene_error_cell = matrix(nrow = 5, ncol = 100)
gene_error_gene = matrix(nrow = 5, ncol = 100)
for(i in c(1:3)){
  error = error_list[[i]]
  m1 = rowMeans(error)
  s1 = apply(error, 1, sd)
  
  error_cell = error_cell_list[[i]]
  m2 = rowMeans(error_cell)
  s2 = apply(error_cell, 1, sd)
  
  error_gene = error_gene_list[[i]]
  m3 = rowMeans(error_gene)
  s3 = apply(error_gene, 1, sd)
  for(j in c(1:5)){
    gene_error[j, ] = rnorm(100, mean = m1[j], sd = s1[j])
    gene_error_cell[j, ] = rnorm(100, mean = m2[j], sd = s2[j])
    gene_error_gene[j, ] = rnorm(100, mean = m3[j], sd = s3[j])
  }
  gene_error_list[[i]] = gene_error
  gene_error_cell_list[[i]] = gene_error_cell
  gene_error_gene_list[[i]] = gene_error_gene
}
saveRDS(gene_error_list, file = "Results/error_all.rds")
saveRDS(gene_error_cell_list, file = "Results/error_all_cell.rds")
saveRDS(gene_error_gene_list, file = "Results/error_all_gene.rds")
