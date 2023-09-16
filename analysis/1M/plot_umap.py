import cudf
#from cuml.manifold.umap import UMAP
import h5py
import pandas as pd
import numpy as np
#dropout
f = h5py.File('/home/suyanchi/project/dab/data/1M/1M.h5')
dat_sc = f['data']
f.close
dat_sc = pd.DataFrame(np.array(dat_sc))

# Initialize UCX for high-speed transport of CUDA arrays
from dask_cuda import LocalCUDACluster

# Create a Dask single-node CUDA cluster w/ one worker per device
cluster = LocalCUDACluster(UDA_VISIBLE_DEVICES="0,1,2,3",
                           protocol="ucx",
                           rmm_pool_size="50GB")
from dask.distributed import Client
client = Client(cluster)
import dask_cudf

from cuml.dask.manifold.umap import UMAP

dat = dask_cudf.from_pandas(dat_sc)

embedding = UMAP().fit_transform(dat)

dask_cudf.Datarame.to_csv('/home/suyanchi/project/dab/results/1M/umap_dropout.csv')
