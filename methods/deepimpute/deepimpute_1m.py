from deepimpute.multinet import MultiNet
import h5py
import pandas as pd
import numpy as np

# load data
f = h5py.File('/home/suyanchi/project/dab/data/1M/1M.h5')
dat_sc = f['data']
f.close
dat_sc = pd.DataFrame(np.array(dat_sc))

model = MultiNet()
model.fit(dat_sc)
imputed = model.predict(dat_sc)

f = h5py.File('/home/suyanchi/project/dab/results/1M/deepimpute.h5', 'w')
f['data'] = imputed.T
f.close