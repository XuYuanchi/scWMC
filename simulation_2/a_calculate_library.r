####### Data Analysis ##############
# Calculate the similarity of two datasets
calculate_similarity <- function(data1,data2){
  
  d = cor(c(data1[lower.tri(data1)]),c(data2[lower.tri(data2)]))
  
  return(d)
  
}
# Calculate the error between true data and imputed data
run_error <- function(dropout_index, seed_value){
  
  # Parameter in the function
  # dropout_index: the index of dropout_mid to control the dropout rate
  # seed_value: the random seed
  
  options( warn = -1 )
  # load the simulationd data
  data_simulation = readRDS(file = paste0('data/simulation_data_dropout_index_',
                                          dropout_index,
                                          '_seed_',
                                          seed_value,
                                          '.rds')
  )
  
  # get the index of the genes with nonzero means
  index = rowMeans(data_simulation$data_dropout) > 0
  
  # obtain the true data
  data_true = data_simulation$data_true
  
  data_true = data_true[index,]
  
  # cell-cell correlation of the true data
  data_true_cell = cor(as.matrix((data_true)), method = "pearson")
  
  # gene-gene correlation of the true data
  data_true_gene = cor(t((data_true)), method = "pearson")
  
  # obtain the dropout data
  data_dropout = data_simulation$data_dropout
  
  data_dropout = data_dropout[index,]
  
  # cell-cell correlation of the dropout data
  data_dropout_cell = cor((data_dropout), method = "pearson")
  
  # gene-gene correlation of the dropout data
  data_dropout_gene = cor(t((data_dropout)), method = "pearson")
  
  # load imputed data from ALRA
  data_ALRA = readRDS( file = paste0("ALRA/",
                                            dropout_index, "_",
                                            seed_value,
                                            ".rds"))
  
  data_ALRA = data_ALRA[index,]
  
  # cell-cell correlation of the ALRA imputed data
  data_ALRA_cell = cor((data_ALRA), method = "pearson")
  
  # gene-gene correlation of the DrImpute imputed data
  data_ALRA_gene = cor(t((data_ALRA)), method = "pearson")
  
  # load imputed data from DCA
  data_DCA = readRDS( file = paste0("DCA/",
                                            dropout_index, "_",
                                            seed_value,
                                            ".rds"))
  
  data_DCA = data_DCA[index,]
  
  # cell-cell correlation of the scImpute imputed data
  data_DCA_cell = cor((data_DCA), method = "pearson")
  
  # gene-gene correlation of the scImpute imputed data
  data_DCA_gene = cor(t((data_DCA)), method = "pearson")
  
  # load the EnImpute imputed data 
  data_EnImpute = readRDS(file = paste0("EnImpute/",dropout_index,"_",seed_value,".rds"))[[1]]
  # browser()
  data_EnImpute[data_EnImpute < 0] = 0

  
  data_EnImpute = data_EnImpute[index,]
  
  # cell-cell correlation of the MAGIC imputed data
  data_EnImpute_cell = cor((data_EnImpute), method = "pearson")
  
  # gene-gene correlation of the MAGIC imputed data
  data_EnImpute_gene = cor(t((data_EnImpute)), method = "pearson")
  
  
  # load imputed data from PBLR
  data_PBLR <- readMat(paste0("PBLR/", dropout_index, "_", seed_value, ".mat"))[[1]]
  
  data_PBLR = data_PBLR[index,]
  
  # cell-cell correlation of the VIPER imputed data
  data_PBLR_cell = cor((data_PBLR), method = "pearson")
  
  # gene-gene correlation of the VIPER imputed data
  data_PBLR_gene = cor(t((data_PBLR)), method = "pearson")
  
  # load imputed data from scVI
  data_scVI = as.matrix(readRDS(file = paste0("scVI/",
                                            dropout_index, "_",
                                            seed_value,
                                            ".rds")))
  data_scVI = data_scVI[index,]

  # cell-cell correlation of the SAVER imputed data
  data_scVI_cell = cor((data_scVI), method = "pearson")
  
  # gene-gene correlation of the SAVER imputed data
  data_scVI_gene = cor(t((data_scVI)), method = "pearson")
  
  # calulate the error between the imputed data and true data
  
  error = matrix(0, nrow = 6, ncol = 1)
  
  error[1] = norm(data_ALRA - log10(data_true + 1), type = "2")
  
  error[2] = norm(log10(data_DCA + 1) - log10(data_true + 1), type = "2")
  
  error[3] = norm(data_EnImpute - log10(data_true + 1), type = "2")
  
  error[4] = norm(data_PBLR - log10(data_true + 1), type = "2")

  error[5] = norm(data_scVI - log10(data_true + 1), type = "2")
  
  error[6] = data_simulation$percentage_zeros
  
  # calulate the similarity between the cell-cell
  # correlation of the imputed data and the one of the true data
  
  error_cell = matrix(0, nrow = 6, ncol = 1)
  
  error_cell[1] = calculate_similarity(data_true_cell, data_ALRA_cell)
  
  error_cell[2] = calculate_similarity(data_true_cell, data_DCA_cell)
  
  error_cell[3] = calculate_similarity(data_true_cell, data_EnImpute_cell)
  
  error_cell[4] = calculate_similarity(data_true_cell, data_PBLR_cell)
  
  error_cell[5] = calculate_similarity(data_true_cell, data_scVI_cell)
  
  error_cell[6] = data_simulation$percentage_zeros
  
  
  # calulate the similarity between the gene-gene
  # correlation of the imputed data and the one of the true da
  error_gene = matrix(0, nrow = 6, ncol = 1)
  
  error_gene[1] = calculate_similarity(data_true_gene, data_ALRA_gene)
  
  error_gene[2] = calculate_similarity(data_true_gene, data_DCA_gene)
  
  error_gene[3] = calculate_similarity(data_true_gene, data_EnImpute_gene)
  
  error_gene[4] = calculate_similarity(data_true_gene, data_PBLR_gene)
  
  error_gene[5] = calculate_similarity(data_true_gene, data_scVI_gene)
  
  error_gene[6] = data_simulation$percentage_zeros
  
  # gather the errors as a list
  result = list()
  
  result$error = error
  
  result$error_cell = error_cell
  
  result$error_gene = error_gene
  
  return(result)
  
}