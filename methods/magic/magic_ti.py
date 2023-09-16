import os
import magic
import pandas as pd
import numpy as np
import h5py
import scipy.io as scio

cwd = os.getcwd()

# load data
# load data
dat_sc = pd.read_csv('/home/suyanchi/project/dab/data/ti/gse.csv', delimiter=',', header=0, index_col=0)
dat_sc = dat_sc.T

magic_operator = magic.MAGIC()
X_magic = magic_operator.fit_transform(dat_sc)


path = ('/home/suyanchi/project/dab/results/ti/magic_gse.mat')
scio.savemat(path, {'re':X_magic.T})