from deepimpute.multinet import MultiNet
import h5py
import pandas as pd
import numpy as np
import scipy.io as scio
import scanpy as sc

dataset_name = ["cell_lines", "panc8_rm", "uc3", "crc", "human_pancreas"]
for i in [4]:
    adata = sc.read_h5ad('/home/suyanchi/project/dab/data/batch/'+dataset_name[i]+'.h5ad')

    model = MultiNet()
    model.fit(adata.to_df())
    # dimension = (cells x genes)
    imputed = model.predict(adata.to_df())
    # path = ('/home/suyanchi/project/dab/results/ti/deepimpute_gse.mat')
    # scio.savemat(path, {'re':imputed.T})
    path = ('/home/suyanchi/project/dab/results/batch/'+ dataset_name[i]+ '_deepimpute.h5')
    f = h5py.File(path, 'w')
    f['data'] = imputed.T
    f.close