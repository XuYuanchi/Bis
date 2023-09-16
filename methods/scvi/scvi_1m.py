import scvi
import scanpy as sc
import numpy as np
import h5py

f = h5py.File('/home/suyanchi/project/dab/data/1M/1M.h5')
dat_sc = f['data']
f.close
dat_sc = np.array(dat_sc)
adata = sc.AnnData(dat_sc)

scvi.model.SCVI.setup_anndata(adata)
model = scvi.model.SCVI(adata)
model.train()
denoised = model.get_normalized_expression(return_numpy=True)

f = h5py.File('/home/suyanchi/project/dab/results/1M/scvi.h5', 'w')
f['data'] = denoised.T
f.close