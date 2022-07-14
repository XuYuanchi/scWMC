color <- readRDS("D:/MyWorkWorld/Imputation/scWMC/color.rds")

get_cor_data <- function(dropout_index, seed_value) {
    options(warn = -1)

    # load the simulationd data

    data_simulation <- readRDS(file = paste0(
        "simulation_data/simulation_data_dropout_index_",
        dropout_index,
        "_seed_",
        seed_value,
        ".rds"
    ))

    index <- rowMeans(data_simulation$data_dropout) > 0

    data_true <- data_simulation$data_true

    data_true <- data_true[index, ]

    # cell-cell correlation

    data_true_cell <- cor(as.matrix((data_true)), method = "pearson")

    data_true_cell[is.na(data_true_cell)] <- 0

    # gene-gene correlation

    data_true_gene <- cor(t((data_true)), method = "pearson")

    data_true_gene[is.na(data_true_gene)] <- 0

    data_dropout <- data_simulation$data_dropout

    data_dropout <- data_dropout[index, ]

    # cell-cell correlation

    data_dropout_cell <- cor((data_dropout), method = "pearson")

    data_dropout_cell[is.na(data_dropout_cell)] <- 0

    # gene-gene correlation

    data_dropout_gene <- cor(t((data_dropout)), method = "pearson")

    data_dropout_gene[is.na(data_dropout_gene)] <- 0


    # load ALRA
    data_ALRA = readRDS(file = paste0(
        "ALRA/",
        dropout_index, "_",
        seed_value,
        ".rds"
    ))

    data_ALRA = data_ALRA[index, ]

    # cell-cell correlation

    data_ALRA_cell <- cor((data_ALRA), method = "pearson")

    data_ALRA_cell[is.na(data_ALRA_cell)] <- 0

    # gene-gene correlation

    data_ALRA_gene <- cor(t((data_ALRA)), method = "pearson")

    data_ALRA_gene[is.na(data_ALRA_gene)] <- 0

    # load DCA
    data_DCA = readRDS(file = paste0(
        "DCA/",
        dropout_index, "_",
        seed_value,
        ".rds"
    ))

    data_DCA = data_DCA[index, ]

    # cell-cell correlation of the scImpute imputed data
    data_DCA_cell = cor((data_DCA), method = "pearson")

    # gene-gene correlation of the scImpute imputed data
    data_DCA_gene = cor(t((data_DCA)), method = "pearson")

    # load the EnImpute imputed data
    data_EnImpute = readRDS(file = paste0("EnImpute/", dropout_index, "_", seed_value, ".rds"))[[1]]
    # browser()
    data_EnImpute[data_EnImpute < 0] = 0


    data_EnImpute = data_EnImpute[index, ]

    # cell-cell correlation of the MAGIC imputed data
    data_EnImpute_cell = cor((data_EnImpute), method = "pearson")

    # gene-gene correlation of the MAGIC imputed data
    data_EnImpute_gene = cor(t((data_EnImpute)), method = "pearson")


    # load imputed data from PBLR
    data_PBLR <- readMat(paste0("PBLR/", dropout_index, "_", seed_value, ".mat"))[[1]]

    data_PBLR = data_PBLR[index, ]

    # cell-cell correlation of the VIPER imputed data
    data_PBLR_cell = cor((data_PBLR), method = "pearson")

    # gene-gene correlation of the VIPER imputed data
    data_PBLR_gene = cor(t((data_PBLR)), method = "pearson")

    # load imputed data from scVI
    data_scVI = as.matrix(readRDS(file = paste0(
        "scVI/",
        dropout_index, "_",
        seed_value,
        ".rds"
    )))
    data_scVI = data_scVI[index, ]

    # cell-cell correlation of the SAVER imputed data
    data_scVI_cell = cor((data_scVI), method = "pearson")

    # gene-gene correlation of the SAVER imputed data
    data_scVI_gene = cor(t((data_scVI)), method = "pearson")
    # load the magic results

    data <- read.table(
        file = paste0("imputation_magic_data/magic_", dropout_index, "_seed_", seed_value, ".csv"),
        header = FALSE,
        sep = ","
    )

    # data$V1 = NULL

    data_magic <- as.matrix(data)

    data_magic[data_magic < 0] <- 0

    data_magic[is.nan(data_magic)] <- 0

    data_magic <- data_magic[index, ]

    # cell-cell correlation

    data_magic_cell <- cor((data_magic), method = "pearson")

    data_magic_cell[is.na(data_magic_cell)] <- 0

    # gene-gene correlation

    data_magic_gene <- cor(t((data_magic)), method = "pearson")

    data_magic_gene[is.na(data_magic_gene)] <- 0

    # load imputed data from scImpute

    data_scimpute <- read.table(
        file = paste0(
            "imputation_scimpute_data/data_imputation_scimpute_",
            dropout_index, "_",
            seed_value,
            ".csv"
        ),
        header = FALSE, sep = ","
    )

    data_scimpute <- data_scimpute[index, ]

    # cell-cell correlation

    data_scimpute_cell <- cor((data_scimpute), method = "pearson")

    data_scimpute_cell[is.na(data_scimpute_cell)] <- 0

    # gene-gene correlation

    data_scimpute_gene <- cor(t((data_scimpute)), method = "pearson")

    data_scimpute_gene[is.na(data_scimpute_gene)] <- 0

    # load imputed data from Drimpute

    data_drimpute <- read.table(
        file = paste0(
            "imputation_drimpute_data/drimpute_",
            dropout_index, "_",
            seed_value,
            ".csv"
        ),
        header = FALSE, sep = ","
    )

    data_drimpute <- data_drimpute[index, ]

    # cell-cell correlation

    data_drimpute_cell <- cor((data_drimpute), method = "pearson")

    data_drimpute_cell[is.na(data_drimpute_cell)] <- 0

    # gene-gene correlation

    data_drimpute_gene <- cor(t((data_drimpute)), method = "pearson")

    data_drimpute_gene[is.na(data_drimpute_gene)] <- 0


    # load imputed data from Saver
    data_viper <- read.table(
        file = paste0(
            "imputation_viper_data/viper_",
            dropout_index, "_",
            seed_value,
            ".csv"
        ),
        header = FALSE, sep = ","
    )

    data_viper <- data_viper[index, ]

    # cell-cell correlation

    data_viper_cell <- cor((data_viper), method = "pearson")

    data_viper_cell[is.na(data_viper_cell)] <- 0

    # gene-gene correlation

    data_viper_gene <- cor(t((data_viper)), method = "pearson")

    data_viper_gene[is.na(data_viper_gene)] <- 0

    # load imputed data from saver

    data_saver <- read.table(
        file = paste0(
            "imputation_saver_data/saver_",
            dropout_index, "_",
            seed_value,
            ".csv"
        ),
        header = FALSE, sep = ","
    )

    data_saver <- data_saver[index, ]


    # cell-cell correlation

    data_saver_cell <- cor(as.matrix(data_saver), method = "pearson")

    data_saver_cell[is.na(data_saver_cell)] <- 0

    # gene-gene correlation

    data_saver_gene <- cor(t((data_saver)), method = "pearson")

    data_saver_gene[is.na(data_saver_gene)] <- 0

    # load imputed data from scTSSR

    data_scTSSR <- read.table(
        file = paste0(
            "imputation_scTSSR_data/scTSSR_",
            dropout_index, "_",
            seed_value,
            ".csv"
        ),
        header = FALSE, sep = ","
    )

    data_scTSSR <- data_scTSSR[index, ]


    # cell-cell correlation

    data_scTSSR_cell <- cor(as.matrix(data_scTSSR), method = "pearson")

    data_scTSSR_cell[is.na(data_scTSSR_cell)] <- 0

    # gene-gene correlation

    data_scTSSR_gene <- cor(t((data_scTSSR)), method = "pearson")

    data_scTSSR_gene[is.na(data_scTSSR_gene)] <- 0
    # load imputed data from our

    data_our <- read.table(
        file = paste0(
            "imputation_our_data/our_",
            dropout_index, "_seed_",
            seed_value,
            ".csv"
        ),
        header = FALSE, sep = ","
    )

    data_our <- data_our[index, ]

    # cell-cell correlation
    data_our_cell <- cor(as.matrix(data_our), method = "pearson")

    data_our_cell[is.na(data_our_cell)] <- 0

    # gene-gene correlation
    data_our_gene <- cor(t((data_our)), method = "pearson")

    data_our_gene[is.na(data_our_gene)] <- 0

    data_cell <- list()

    data_cell[[1]] <- data_true_cell

    data_cell[[2]] <- data_dropout_cell

    data_cell[[3]] <- data_our_cell

    data_cell[[4]] <- data_ALRA_cell

    data_cell[[5]] <- data_DCA_cell

    data_cell[[6]] <- data_drimpute_gene

    data_cell[[7]] <- data_EnImpute_cell

    data_cell[[8]] <- data_magic_cell

    data_cell[[9]] <- data_PBLR_cell

    data_cell[[10]] <- data_saver_cell

    data_cell[[11]] <- data_scimpute_cell

    data_cell[[12]] <- data_scTSSR_cell

    data_cell[[13]] <- data_scVI_cell

    data_cell[[14]] <- data_viper_cell


    data_gene <- list()

    data_gene[[1]] <- data_true_gene

    data_gene[[2]] <- data_dropout_gene

    data_gene[[3]] <- data_our_gene

    data_gene[[4]] <- data_ALRA_gene

    data_gene[[5]] <- data_DCA_gene

    data_gene[[6]] <- data_drimpute_gene

    data_gene[[7]] <- data_EnImpute_gene

    data_gene[[8]] <- data_magic_gene

    data_gene[[9]] <- data_PBLR_gene

    data_gene[[10]] <- data_saver_gene

    data_gene[[11]] <- data_scimpute_gene

    data_gene[[12]] <- data_scTSSR_gene

    data_gene[[13]] <- data_scVI_gene

    data_gene[[14]] <- data_viper_gene

    data_cor <- list()

    data_cor[[1]] <- data_cell

    data_cor[[2]] <- data_gene

    return(data_cor)
}

data_cor_vector <- function(data) {
    return(data[lower.tri(data)])
}

plot_cell_distribution <- function(dropout_index, seed_value) {
    methods <- c(
        "True Data", "Dropout Data", "scWMC", "ALRA", "DCA", "DrImpute", "EnImpute", 
        "MAGIC", "PBLR", "SAVER", "scImpute", "scTSSR", "scVI", "VIPER"
    )
    # load the simulationd data
    data_simulation <- readRDS(file = paste0(
        "simulation_data/simulation_data_dropout_index_",
        dropout_index,
        "_seed_",
        seed_value,
        ".rds"
    ))

    group_info <- unique(data_simulation$group)
    #browser()

    # data_cor <- get_cor_data(dropout_index, seed_value)

    # saveRDS(data_cor, paste0("D:/MyWorkWorld/Imputation/scWMC/simulation_1/cor/", dropout_index, "_", seed_value, ".rds"))
    
    data_cor = readRDS("D:/MyWorkWorld/Imputation/scWMC/simulation_1/cor/1_10.rds")

    data_cell <- data_cor[[1]]

    # data_gene <- data_cor[[2]]

    k <- 1

    p <- list()

    for (j in group_info) {
        index1 <- data_simulation$group == j

        ratio <- c()

        for (i in c(1:length(methods))) {
            tmp <- data_cell[[i]]
            
            index1 = index1[1:dim(tmp)[1]]

            low_index <- lower.tri(tmp) * 1

            tmp_index <- matrix(0, nrow = dim(tmp)[1], ncol = dim(tmp)[2])

            tmp_index[index1, index1] <- 1

            low_index1 <- low_index * tmp_index > 0

            tmp_index <- matrix(0, nrow = dim(tmp)[1], ncol = dim(tmp)[2])

            tmp_index[index1, !index1] <- 1

            tmp_index[!index1, index1] <- 1

            low_index2 <- low_index * tmp_index > 0

            data_plot1 <- data.frame(x = data_cor_vector(tmp[low_index1]))

            data_plot1$y <- j

            data_plot2 <- data.frame(x = data_cor_vector(tmp[low_index2]))

            data_plot2$y <- "Group0"

            data_plot <- rbind(data_plot1, data_plot2)

            a <- ks.test(data_plot1$x, data_plot2$x)

            ratio <- cbind(ratio, a[[1]])

            p[[k]] <- ggplot(data_plot, aes(x = x, fill = y)) +
                geom_density(alpha = 0.4) +
                ggtitle(paste0("Cell: ", methods[i])) +
                theme_bw() +
                theme(
                  axis.title.x=element_blank(),
                    legend.position = "bottom",
                    plot.title = element_text(face = "bold", size = 14, hjust = 0.4),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank()
                ) +
                scale_fill_manual(values = c("#DC0000FF", "#3C5488FF"))

            k <- k + 1
        }
        methods <- factor(methods, levels = c(
            "True Data", "Dropout Data", "scWMC", "ALRA", "DCA", "DrImpute", "EnImpute",
            "MAGIC", "PBLR", "SAVER", "scImpute", "scTSSR", "scVI", "VIPER"
        ))
        data_bar <- data.frame(y = t(ratio), x = methods)

        colnames(data_bar) <- c("y", "x")

        p[[k]] <- ggplot(data_bar, aes(x = x, y = y, fill = methods)) +
            geom_bar(stat = "identity") +
            # scale_x_discrete(limits=data_bar$x) +
            xlab(NULL) +
            ylab("KS Statistics") +
            theme_bw() +
            scale_fill_manual(values = c("#7E6148FF", color)) +
            scale_x_discrete(labels = c("True Data", "Dropout Data", "scWMC", "ALRA", "DCA", "DrImpute", "EnImpute",
            "MAGIC", "PBLR", "SAVER", "scImpute", "scTSSR", "scVI", "VIPER")) +
            # ggtitle(paste0(j, "_", round(data_simulation$percentage_zeros))) +
            # theme_cowplot() +
            theme(
                plot.title = element_text(face = "bold", size = 14, hjust = 0.4),
                axis.text.x = element_text(angle = -30, size = 10, face = "bold"),
                legend.position = "bottom",
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank()
            )



        k <- k + 1
    }

    # main = grid.arrange(grobs = p,ncol = (length(methods) + 1))

    return(p)
}

plot_gene_distribution <- function(dropout_index, seed_value) {
    methods <- c(
        "True Data", "Dropout Data", "scWMC", "ALRA", "DCA", "DrImpute", "EnImpute",
        "MAGIC", "PBLR", "SAVER", "scImpute", "scTSSR", "scVI", "VIPER"
    )
    # load the simulationd data
    data_simulation <- readRDS(file = paste0(
        "simulation_data/simulation_data_dropout_index_",
        dropout_index,
        "_seed_",
        seed_value,
        ".rds"
    ))

    group_info <- unique(data_simulation$group)

    data_cor <- readRDS(paste0("D:/MyWorkWorld/Imputation/scWMC/simulation_1/cor/", dropout_index, "_", seed_value, ".rds"))

    # data_cell <- data_cor[[1]]

    data_gene <- data_cor[[2]]

    tmp <- as.numeric(as.numeric(gsub("Group", "", data_simulation$group)))

    de <- get_marker_genes(data_simulation$data_true, tmp)

    group_unique <- unique(tmp)

    N_groups <- length(group_unique)

    de_gene <- list()

    for (i in c(1:N_groups)) {
        de_gene[[i]] <- paste0("Gene", which((de$auroc > 0.85) & (de$clusts == group_unique[i]) & (de$pvalue < 0.01)))
    }

    names_gene <- rownames(data_gene[[1]])

    k <- 1

    p <- list()

    for (j in c(1:N_groups)) {
        index1 <- names_gene %in% de_gene[[j]]

        ratio <- c()

        for (i in c(1:length(methods))) {
            tmp <- data_gene[[i]]

            low_index <- lower.tri(tmp) * 1

            tmp_index <- matrix(0, nrow = dim(tmp)[1], ncol = dim(tmp)[2])

            tmp_index[index1, index1] <- 1

            low_index1 <- low_index * tmp_index > 0

            tmp_index <- matrix(0, nrow = dim(tmp)[1], ncol = dim(tmp)[2])

            tmp_index[index1, !index1] <- 1

            tmp_index[!index1, index1] <- 1

            low_index2 <- low_index * tmp_index > 0

            data_plot1 <- data.frame(x = data_cor_vector(tmp[low_index1]))

            data_plot1$y <- paste0("Marker_", j)

            data_plot2 <- data.frame(x = data_cor_vector(tmp[low_index2]))

            data_plot2$y <- "NonMarker"

            data_plot <- rbind(data_plot1, data_plot2)

            a <- ks.test(data_plot1$x, data_plot2$x)

            ratio <- cbind(ratio, a[[1]])

            p[[k]] <- ggplot(data_plot, aes(x = x, fill = y)) +
                geom_density(alpha = 0.7) +
                ggtitle(paste0("Gene: ", methods[i])) +
                theme_bw() +
                theme(
                    axis.title.x=element_blank(),
                    plot.title = element_text(face = "bold", size = 14, hjust = 0.4),
                    legend.position = "bottom",
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank()
                ) +
                scale_fill_manual(values = c("#DC0000FF", "#3C5488FF"))

            k <- k + 1
        }

        methods <- factor(methods, levels = c(
            "True Data", "Dropout Data", "scWMC", "ALRA", "DCA", "DrImpute", "EnImpute",
            "MAGIC", "PBLR", "SAVER", "scImpute", "scTSSR", "scVI", "VIPER"
        ))
        data_bar <- data.frame(y = t(ratio), x = methods)

        colnames(data_bar) <- c("y", "x")

        p[[k]] <- ggplot(data_bar, aes(x = x, y = y, fill = methods)) +
            geom_bar(stat = "identity") +
            # scale_x_discrete(limits=data_bar$x) +
            xlab(NULL) +
            ylab("KS Statistics") +
            # ggtitle(paste0(j, "_", round(data_simulation$percentage_zeros))) +
            # theme_cowplot() +
            scale_fill_manual(values = c("#7E6148FF", color)) +
            scale_x_discrete(labels = c("True Data", "Dropout Data", "scWMC", "ALRA", "DCA", "DrImpute", "EnImpute",
        "MAGIC", "PBLR", "SAVER", "scImpute", "scTSSR", "scVI", "VIPER")) +
            theme_bw() +
            theme(
                plot.title = element_text(face = "bold", size = 14, hjust = 0.4),
                axis.text.x = element_text(angle = -30, size = 10, face = "bold"),
                legend.position = "bottom",
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank()
            )

        k <- k + 1
    }

    # main = grid.arrange(grobs = p,ncol = (length(methods) + 1))

    return(p)
}