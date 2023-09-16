from dca.api import dca
import h5py
import numpy as np
import pandas as pd
import scanpy as sc
import time

# f = h5py.File('/home/suyanchi/project/dab/data/1M/1M.h5')
f = h5py.File('/home/suyanchi/project/dab/data/time/1000.h5')
dat_sc = f['data']
f.close
t1 = time.time()
dat_sc = pd.DataFrame(np.array(dat_sc)[:,1:1000])

adata = sc.AnnData(dat_sc.values)
print(adata.shape)
result = dca(adata, copy=True)

t2 = time.time()
print(t2-t1)