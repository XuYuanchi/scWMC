import scvi
import scanpy as sc
import numpy as np
import torch
data_name = ["sc_CELseq2", "sc_10x", "sc_Droseq", "sc_10x_5cl", "sc_Celseq2_5cl_p1", "sc_Celseq2_5cl_p2", "sc_Celseq2_5cl_p3"]
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
for j in range(1,11):
    for i in range(2,6):
        filename = "/home/suyanchi/project/scWMC/Clustering/data/" + data_name[i] + ".csv"
        adata = scvi.data.read_csv(filename, delimiter=',', first_column_names=True, dtype='float32')
        adata = adata.transpose()
        scvi.model.SCVI.setup_anndata(adata)
        model = scvi.model.SCVI(adata)
        model.to_device(device)
        model.train()
        denoised = model.get_normalized_expression(return_numpy=True)
        np.savetxt("/home/suyanchi/project/scWMC/Clustering/scVI/" + str(j) + "/" + data_name[i] + ".csv", denoised, delimiter=",")
