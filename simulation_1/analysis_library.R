######## Data Generation ##############
# the function generating the simulation data using bioconductor package Splatter
generate_simulation_splatter <- function(dropout_index, seed_value, nGenes = 800){
  
  # Parameter in the function
  # dropout_index: the index of dropout_mid to control the dropout rate
  #    seed_value: the random seed
  #        nGenes: the number of genes in the simulation data. The default is 800
  
  # Set up the parameters
  params = newSplatParams()
  
  params = setParams(params, list(batchCells = 1000,
                                  nGenes = nGenes,
                                  group.prob = c(0.20, 0.35, 0.45),
                                  de.prob = c(0.045, 0.045, 0.045),
                                  de.facLoc = 0.1,
                                  de.facScale = 0.4)
  )
  
  # Set up the vector of dropout.mid
  dropout_mid = c(4, 5, 5.5)
  
  # determine if it is a good parameter
  if(dropout_index > length(dropout_mid)){
    
    stop(
      paste0('The dropout_index shold not be greater than ', 
             length(dropout_mid), 
             ' . Please input a proper one.\n')
    )
    
  }
  
  # Generate the simulation data using Splatter package
  sim = splatSimulateGroups(params,
                            dropout.type = "experiment",
                            dropout.shape = -1,
                            dropout.mid = dropout_mid[dropout_index],
                            seed = seed_value)
  # browser() 
  # genereate the cpm levels of the true simulation data
  data_true = cpm(sim@assays@data$TrueCounts)
  
  data_dropout = data_true
  browser()
  a = (counts(sim) == 0)
  if(class(a)=='matrix'){
    data_dropout[a] = 0
    }
  else{
    data_dropout[a@x] = 0
  }
  # generate the dropout data based on the counts in sim
  # data_dropout[(counts(sim) == 0)] = 0
  
  # calculate the dropout rate
  percentage_zeros = round(nnzero(data_dropout == 0, na.counted = NA)/
                             (dim(data_dropout)[1]*dim(data_dropout)[2])*100)
  
  
  # generate the bulk RNAseq data
  data_bulk = data.frame(val = rowMeans(data_true))
  
  # define the data list for the simulation data
  # indcluding: data_true: true data
  #          data_dropout: dropout data
  #             data_bluk: bulk data
  #      percentage_zeros: dropout rate
  #                 group: the group label
  
  data = list()
  
  data$data_bulk = data_bulk
  
  data$data_dropout = data_dropout
  
  data$data_true = data_true
  
  data$percentage_zeros = percentage_zeros
  
  data$group = colData(sim)@listData$Group
  
  return(data)
}

# generate the simulation data and save the data
generate_save_data <- function(dropout_index, seed_value) {

  # Parameter in the function
  # dropout_index: the index of dropout_mid to control the dropout rate
  # seed_value: the random seed

  # generate the simulation data
  data_simulation = generate_simulation_splatter(dropout_index, seed_value)

  # generate the folder saving the simulation data
  dir.create(file.path("simulation_data"), showWarnings = FALSE)

  # save the data as RDS format
  saveRDS(data_simulation,
    file = paste0(
      "simulation_data/simulation_data_dropout_index_",
      dropout_index,
      "_seed_",
      seed_value,
      ".rds"
    )
  )
  filename <- paste("simulation_data/simulation_data_dropout_index_", dropout_index,
    "_seed_", seed_value, ".mat",
    sep = ""
  )
  writeMat(filename,
    data_bulk = data_simulation$data_bulk, data_dropout = data_simulation$data_dropout,
    data_true = data_simulation$data_true
  )
}