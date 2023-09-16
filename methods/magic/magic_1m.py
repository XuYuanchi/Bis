import os
import magic
import pandas as pd
import numpy as np
import h5py

cwd = os.getcwd()

# load data
f = h5py.File('/home/suyanchi/project/dab/data/1M/1M.h5')
dat_sc = f['data']
f.close
dat_sc = pd.DataFrame(np.array(dat_sc))

magic_operator = magic.MAGIC()
X_magic = magic_operator.fit_transform(dat_sc)
out_magic = X_magic.T

f = h5py.File('/home/suyanchi/project/dab/results/1M/magic.h5', 'w')
f['data'] = out_magic
f.close