import scvi
import scanpy as sc
import numpy as np
import torch
data_name = ["Deng", "Petropoulos"]
device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
  
filename = 'D:/MyWorkWorld/Imputation/scWMC/ICT/data/' + data_name[1] + '_filter.csv'
adata = scvi.data.read_csv(filename, delimiter=',', first_column_names=True, dtype='float32')
adata = adata.transpose()
scvi.model.SCVI.setup_anndata(adata)
model = scvi.model.SCVI(adata)
model.to_device(device)
model.train()
denoised = model.get_normalized_expression(return_numpy=True)
np.savetxt('D:/MyWorkWorld/Imputation/scWMC/ICT/Results/' + data_name[1] + '/scVI.csv', denoised, delimiter=',')
