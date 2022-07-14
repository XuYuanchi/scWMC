import scvi
import scanpy as sc
import numpy as np

filename = "/home/grads/fuzhowang2/suyanchi/scWMC/scalability/data/100k.csv"
adata = scvi.data.read_csv(filename, delimiter=',', first_column_names=True, dtype='float32')
scvi.model.SCVI.setup_anndata(adata)
adata = adata.transpose()
model = scvi.model.SCVI(adata)
model.train()
denoised = model.get_normalized_expression(return_numpy=True)
np.savetxt("/home/grads/fuzhowang2/suyanchi/scWMC/scalability/Results/scVI.csv", denoised, delimiter=",")