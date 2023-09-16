import os
import magic
import pandas as pd
import numpy as np
import h5py
import scipy.io as scio
import scanpy as sc

dataset_name = ["cell_lines", "panc8_rm", "uc3", "crc", "human_pancreas"]
for i in [4]:
    adata = sc.read_h5ad('/home/suyanchi/project/dab/data/batch/'+dataset_name[i]+'.h5ad')

    magic_operator = magic.MAGIC()
    # dimension = (cells x genes)
    X_magic = magic_operator.fit_transform(adata.to_df().values)

    path = '/home/suyanchi/project/dab/results/batch/'+ dataset_name[i]+ '_magic.h5'
    f = h5py.File(path, 'w')
    f['data'] = X_magic
    f.close