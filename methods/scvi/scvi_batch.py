import scvi
import scanpy as sc
import numpy as np
import pandas as pd
import h5py
import scipy.io as scio
import torch

import os
os.environ["CUDA_VISIBLE_DEVICES"] = "2"


device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

dataset_name = ["cell_lines", "panc8_rm", "uc3", "crc", "human_pancreas"]
for i in [4]:
    adata = sc.read_h5ad('/home/suyanchi/project/dab/data/batch/'+dataset_name[i]+'.h5ad')

    scvi.model.SCVI.setup_anndata(adata)
    model = scvi.model.SCVI(adata)
    model.to_device(device)
    model.train()
    denoised = model.get_normalized_expression(return_numpy=True)

    path = ('/home/suyanchi/project/dab/results/batch/'+ dataset_name[i]+ '_scvi.h5')
    f = h5py.File(path, 'w')
    f['data'] = denoised
    f.close