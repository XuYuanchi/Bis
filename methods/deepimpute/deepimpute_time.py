from deepimpute.multinet import MultiNet
import h5py
import pandas as pd
import numpy as np
import time

# load data
# f = h5py.File('/home/suyanchi/project/dab/data/time/1000000.h5')
f = h5py.File('/home/suyanchi/project/dab/data/1M/1M.h5')
dat_sc = f['data']
f.close
t1 = time.time()
dat_sc = pd.DataFrame(np.array(dat_sc)[:,1:1000])

model = MultiNet()
model.fit(dat_sc)
imputed = model.predict(dat_sc)
t2 = time.time()
print(t2-t1)
