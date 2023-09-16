from dca.api import dca
import h5py
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io as scio

dataset_name = ["cell_lines", "panc8_rm", "uc3", "crc", "human_pancreas"]
for i in [4]:
    adata = sc.read_h5ad('/home/suyanchi/project/dab/data/batch/'+dataset_name[i]+'.h5ad')

    result = dca(adata, copy=True, check_counts=False)

    path = '/home/suyanchi/project/dab/results/batch/'+ dataset_name[i]+ '_dca.h5'
    f = h5py.File(path, 'w')
    f['data'] = result.X
    f.close
