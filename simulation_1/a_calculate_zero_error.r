library(R.matlab)

for (dropout_index in c(1:3)) {
    for (seed_value in c(1:10)) {

        # read raw
        data_simulation = readRDS(file = paste0(
            "D:/project/Demo_4/SimulationData_Strategy1/simulation_data/simulation_data_dropout_index_",
            dropout_index, "_seed_", seed_value, ".rds"
        ))
        # get the index of the genes with nonzero means
        index = rowMeans(data_simulation$data_dropout) > 0
        
        data_true    = data_simulation$data_true
        
        data_true = data_true[index,]
        
        data_dropout = data_simulation$data_dropout
        
        data_dropout = data_dropout[index,]
        
        mask         = data_dropout
        mask[mask > 0] = -1
        mask[mask == 0] = 1
        mask[mask == -1] = 0
        data_true = data_true * mask

        # load imputed data from our
        data_our = read.table(
            file = paste0(
                "D:/project/Demo_4/SimulationData_Strategy1/imputation_our_data/our_",
                dropout_index, "_seed_",
                seed_value,
                ".csv"
            ),
            header = FALSE, sep = ","
        )
        data_our = data_our[index,]
        # # load imputed data from ALRA
        data_ALRA = readRDS(file = paste0(
            "D:/project/Demo_4/SimulationData_Strategy1/ALRA/",
            dropout_index, "_",
            seed_value,
            ".rds"
        ))
        data_ALRA = data_ALRA[index,]
        # load imputed data from DCA
        data_DCA = readRDS(file = paste0(
            "D:/project/Demo_4/SimulationData_Strategy1/DCA/",
            dropout_index, "_",
            seed_value,
            ".rds"
        ))
        data_DCA = data_DCA[index,]
        # load imputed data from DrImpute
        data_drimpute = read.table(
            file = paste0(
                "D:/project/Demo_4/SimulationData_Strategy1/imputation_drimpute_data/drimpute_",
                dropout_index, "_",
                seed_value,
                ".csv"
            ),
            header = FALSE, sep = ","
        )
        data_drimpute = data_drimpute[index,]
        # load the EnImpute imputed data
        data_EnImpute = readRDS(file = paste0("D:/project/Demo_4/SimulationData_Strategy1/EnImpute/", dropout_index, "_", seed_value, ".rds"))[[1]]
        data_EnImpute = data_EnImpute[index,]
        # load the MAGIC imputed data
        data_magic = read.table(
            file = paste0("D:/project/Demo_4/SimulationData_Strategy1/imputation_magic_data/magic_", dropout_index, "_seed_", seed_value, ".csv"),
            header = FALSE,
            sep = ","
        )
        data_magic = data_magic[index,]
        # load imputed data from PBLR
        data_PBLR <- readMat(paste0("D:/project/Demo_4/SimulationData_Strategy1/PBLR/", dropout_index, "_", seed_value, ".mat"))[[1]]
        data_PBLR = data_PBLR[index,]
        # load imputed data from SAVER
        data_saver = read.table(
            file = paste0(
                "D:/project/Demo_4/SimulationData_Strategy1/imputation_saver_data/saver_",
                dropout_index, "_",
                seed_value,
                ".csv"
            ),
            header = FALSE, sep = ","
        )
        data_saver = data_saver[index,]
        # load imputed data from scImpute
        data_scimpute = read.table(
            file = paste0(
                "D:/project/Demo_4/SimulationData_Strategy1/imputation_scimpute_data/data_imputation_scimpute_",
                dropout_index, "_",
                seed_value,
                ".csv"
            ),
            header = FALSE, sep = ","
        )
        data_scimpute = data_scimpute[index,]
         # load imputed data from scTSSR
         data_scTSSR = read.table(
             file = paste0(
                 "D:/project/Demo_4/SimulationData_Strategy1/imputation_scTSSR_data/scTSSR_",
                 dropout_index, "_",
                 seed_value,
                 ".csv"
             ),
             header = FALSE, sep = ","
         )
         data_scTSSR = data_scTSSR[index,]
          # load imputed data from scVI
          data_scVI = as.matrix(readRDS(file = paste0(
              "D:/project/Demo_4/SimulationData_Strategy1/scVI/",
              dropout_index, "_",
              seed_value,
              ".rds"
          )))
          data_scVI = data_scVI[index,]
          # load imputed data from VIPER
          data_viper = read.table(
              file = paste0(
                  "D:/project/Demo_4/SimulationData_Strategy1/imputation_viper_data/viper_",
                  dropout_index, "_",
                  seed_value,
                  ".csv"
              ),
              header = FALSE, sep = ","
          )
          data_viper = data_viper[index,]
          error = matrix(0, nrow = 14, ncol = 1)

          error[1] = log2(norm(data_true, type = "2"))

          error[2] = log2(norm(data_our * mask - data_true, type = "2"))

          error[3] = log2(norm((exp(data_ALRA) - 1) * mask - data_true, type = "2"))

          error[4] = log2(norm(data_DCA * mask - data_true, type = "2"))

          error[5] = log2(norm(data_drimpute * mask - data_true, type = "2"))

          error[6] = log2(norm(data_EnImpute * mask - data_true, type = "2"))

          error[7] = log2(norm(data_magic * mask - data_true, type = "2"))

          error[8] = log2(norm((10^data_PBLR - 1) * mask - data_true, type = "2"))

          error[9] = log2(norm(data_saver * mask - data_true, type = "2"))

          error[10] = log2(norm(data_scimpute * mask - data_true, type = "2"))

          error[11] = log2(norm(data_scTSSR * mask - data_true, type = "2"))

          error[12] = log2(norm(data_scVI * mask - data_true, type = "2"))

          error[13] = log2(norm(data_viper * mask - data_true, type = "2"))

          error[14] = data_simulation$percentage_zeros

          saveRDS(error,
              file = paste0("D:/MyWorkWorld/Imputation/scWMC/simulation_1/error_data/error_zero", dropout_index, "_", seed_value, ".rds"))
          
          
          
          
    }
}

# gather error

error_list <- list()

for (i in c(1:3)) {
    error_matrix <- c()

    for (j in c(1:10)) {
        tmp <- readRDS(file = paste0("D:/MyWorkWorld/Imputation/scWMC/simulation_1/error_data/error_zero", i, "_", j, ".rds"))
        error_matrix <- cbind(error_matrix, as.matrix(tmp))
    }

    error_list[[i]] <- error_matrix
}

saveRDS(error_list, file = "D:/MyWorkWorld/Imputation/scWMC/simulation_1/data/error_zero_all.rds")
