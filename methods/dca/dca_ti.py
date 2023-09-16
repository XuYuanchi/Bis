from dca.api import dca
import h5py
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io as scio

dat_sc = pd.read_csv('/home/suyanchi/project/dab/data/ti/gse.csv', delimiter=',', header=0, index_col=0)

adata = sc.AnnData(dat_sc.T.astype(np.int))
result = dca(adata, copy=True)

path = ('/home/suyanchi/project/dab/results/ti/dca_gse.mat')
scio.savemat(path, {'re':result.X.T})
