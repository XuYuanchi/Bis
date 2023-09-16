from dca.api import dca
import h5py
import numpy as np
import pandas as pd
import scanpy as sc

f = h5py.File('/home/suyanchi/project/dab/data/1M/1M.h5')
dat_sc = f['data']
f.close
dat_sc = pd.DataFrame(np.array(dat_sc))

adata = sc.AnnData(dat_sc.values)
result = dca(adata, copy=True)

f = h5py.File('/home/suyanchi/project/dab/results/1M/dca.h5', 'w')
f['data'] = result.X.T
f.close