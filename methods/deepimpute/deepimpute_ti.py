from deepimpute.multinet import MultiNet
import h5py
import pandas as pd
import numpy as np
import scipy.io as scio

# load data
dat_sc = pd.read_csv('/home/suyanchi/project/dab/data/ti/gse.csv', delimiter=',', header=0, index_col=0)
dat_sc = dat_sc.T

model = MultiNet()
model.fit(dat_sc)
imputed = model.predict(dat_sc)

# path = ('/home/suyanchi/project/dab/results/ti/deepimpute_gse.mat')
# scio.savemat(path, {'re':imputed.T})
f = h5py.File('/home/suyanchi/project/dab/results/ti/deepimpute_gse.h5', 'w')
f['data'] = imputed.T
f.close