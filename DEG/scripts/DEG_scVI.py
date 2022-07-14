import scvi
import scanpy as sc
import numpy as np


filename = "/home/grads/fuzhowang2/suyanchi/scWMC/DEG/data/DEG.raw.csv"
adata = scvi.data.read_csv(filename, delimiter=',', first_column_names=False, dtype='float32')
adata = adata.transpose()
scvi.model.SCVI.setup_anndata(adata)
model = scvi.model.SCVI(adata)
model.train()
denoised = model.get_normalized_expression(return_numpy=True)
np.savetxt("/home/grads/fuzhowang2/suyanchi/scWMC/DEG/Results/scVI.csv", denoised, delimiter=",")
