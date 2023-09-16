import scvi
import scanpy as sc
import numpy as np
import h5py
import time

f = h5py.File('/home/suyanchi/project/dab/data/time/500000.h5')
# f = h5py.File('/home/suyanchi/project/dab/data/1M/1M.h5')
dat_sc = f['data']
f.close
t1 = time.time()
dat_sc = np.array(dat_sc)
adata = sc.AnnData(dat_sc.T)
print(adata.shape)

scvi.model.SCVI.setup_anndata(adata)
model = scvi.model.SCVI(adata)
model.to_device(3)
model.train()
denoised = model.get_normalized_expression(return_numpy=True)
t2 = time.time()
print(t2-t1)

