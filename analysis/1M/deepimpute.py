import h5py
import umap
import numpy as np
import pandas as pd
# from sklearn.preprocessing import StandardScaler

reducer = umap.UMAP()

#dropout
f = h5py.File('/home/suyanchi/project/dab/results/1M/deepimpute.h5')
dat_sc = f['data']
f.close
dat_sc = np.array(dat_sc)
print(dat_sc.shape)
#scaled_dat_sc = StandardScaler().fit_transform(dat_sc)
embedding = reducer.fit_transform(dat_sc.T)

f = h5py.File('/home/suyanchi/project/dab/results/1M/umap_deepimpute.h5', 'w')
f['data'] = embedding
f.close





