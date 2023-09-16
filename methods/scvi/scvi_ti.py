import scvi
import scanpy as sc
import numpy as np
import pandas as pd
import h5py
import scipy.io as scio
import torch

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

dat_sc = pd.read_csv('/home/suyanchi/project/dab/data/ti/gse.csv', delimiter=',', header=0, index_col=0)
adata = sc.AnnData(dat_sc.T)

scvi.model.SCVI.setup_anndata(adata)
model = scvi.model.SCVI(adata)
model.to_device(device)
model.train()
denoised = model.get_normalized_expression(return_numpy=True)

path = ('/home/suyanchi/project/dab/results/ti/scvi_gse.mat')
scio.savemat(path, {'re':denoised.T})