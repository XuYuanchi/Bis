import os
import anndata as ad
import numpy as np
import pandas as pd
from pandas import DataFrame, Series
import pytorch_lightning as pl
import scanpy as sc
from torch.utils.data import Dataset, DataLoader, random_split
import h5py
import scipy.io as sio
import scipy.sparse as sp

class DatasetWithConfounder(Dataset):
    def __init__(self, X: DataFrame, sf: DataFrame):
        self.X = X
        self.sf = sf

    def __len__(self):
        return self.X.shape[0]

    def __getitem__(self, index):
        if type(self.X) == DataFrame:
            X_ret = self.X.iloc[index].to_numpy()
            sf = self.sf.iloc[index].to_numpy()
        else:
            X_ret = self.X[index]
            sf = self.sf[index]
        return [X_ret, sf]


class DataModule(pl.LightningDataModule):
    def __init__(self, data_dir='/home/suyanchi/project/dab/data/', data_name='', batchsize=512):
        super().__init__()
        self.data_dir = data_dir
        self.batch_size = batchsize
        self.ann_data = None
        self.Dataset = data_name

    def prepare_data(self, filter_min_counts=False, size_factors=False, normalize_input=False, logtrans_input=True):
        # loading and proprecessing data
        if self.Dataset == '1M':
            f = h5py.File('/home/suyanchi/project/dab/data/1M/1M.h5')
            dat_sc = f['data']
            f.close
            self.ann_data = sc.AnnData(pd.DataFrame(np.array(dat_sc)).values)
        elif self.Dataset in ["cell_lines", "panc8_rm", "uc3", "crc",'human_pancreas']:
            self.ann_data = sc.read_h5ad('/home/suyanchi/project/dab/data/batch/'+ self.Dataset+'.h5ad')
            if sp.issparse(self.ann_data.X):
                self.ann_data.X = self.ann_data.X.toarray()
            # normalize_input = True
        elif self.Dataset == 'PBMC':
            self.ann_data = sc.read_h5ad('/home/suyanchi/project/dab/data/case/'+ self.Dataset+'.h5ad')
            if sp.issparse(self.ann_data.X):
                self.ann_data.X = self.ann_data.X.toarray()
        elif self.Dataset == 'time':
            # f = h5py.File('/home/suyanchi/project/dab/data/time/500000.h5')
            f = h5py.File('/home/suyanchi/project/dab/data/1M/1M.h5')
            dat_sc = f['data']
            f.close
            self.ann_data = sc.AnnData(pd.DataFrame(np.array(dat_sc[:,1:1000])).values)
            # self.ann_data = sc.AnnData(pd.DataFrame(np.array(dat_sc)).values.T)
        elif self.Dataset == 'test':
            file_path = self.data_dir + self.Dataset + '/6/1.mat'
            f = sio.loadmat(file_path)
            self.ann_data = sc.AnnData(pd.DataFrame(np.array(f['data_dropout']).T))
        elif self.Dataset == 'deg':
            file_path = self.data_dir + self.Dataset + '/sc_2000.csv'
            self.ann_data = sc.read_csv(file_path, first_column_names=True)
            self.ann_data = self.ann_data.T
        elif self.Dataset == 'deg_raw':
            file_path = self.data_dir + self.Dataset + '/sc_raw.csv'
            self.ann_data = sc.read_csv(file_path, first_column_names=True)
            self.ann_data = self.ann_data.T
        elif self.Dataset in ['liver', 'heart', 'marrow', 'lung']:
            file_path = self.data_dir + 'downsample/' + self.Dataset + '.mat'
            f= sio.loadmat(file_path)
            self.ann_data = sc.AnnData(pd.DataFrame(np.array(f['data_sc']).T))
        elif self.Dataset in ['gse', 'Deng', 'Petropoulos']:
            file_path = self.data_dir + 'ti/' + self.Dataset + '.csv'
            self.ann_data = sc.read_csv(file_path, first_column_names=True)
            self.ann_data = self.ann_data.T
        else:
            file_path = self.data_dir + self.Dataset + '_dropout.csv'
            self.ann_data = sc.read_csv(file_path, first_column_names=True)
            self.ann_data = self.ann_data.T
        # filter cell and gene
        if filter_min_counts:
            sc.pp.filter_genes(self.ann_data, min_counts=1)
            sc.pp.filter_cells(self.ann_data, min_counts=1)

        if size_factors or normalize_input or logtrans_input:
            self.ann_data.raw = self.ann_data.copy()
        else:
            self.ann_data.raw = self.ann_data

        if size_factors:
            sc.pp.normalize_per_cell(self.ann_data)
            self.ann_data.obs['size_factors'] = self.ann_data.obs.n_counts / np.median(self.ann_data.obs.n_counts)
        else:
            self.ann_data.obs['size_factors'] = 1.0

        if logtrans_input:
            sc.pp.log1p(self.ann_data)

        if normalize_input:
            sc.pp.scale(self.ann_data)

    def setup(self, stage=None):
        # split datasets
        self.data_set = DatasetWithConfounder(X=self.ann_data.X.copy(), sf=self.ann_data.obs['size_factors'])
        if stage == 'fit' or stage is None:
            size = len(self.data_set)
            t, v = (int(size * 0.9), int(size * 0.1))
            t += (t + v != size)
            self.dataset_train, self.dataset_val = random_split(self.data_set, [t, v])

        if stage == 'test' or stage is None:
            self.dataset_test = self.data_set

    def train_dataloader(self):
        return DataLoader(
            self.dataset_train,
            batch_size=self.batch_size,
            num_workers=3,
            shuffle=False,
            drop_last=False
        )

    def val_dataloader(self):
        return DataLoader(
            self.dataset_val,
            batch_size=self.batch_size,
            num_workers=3
        )

    def test_dataloader(self):
        return DataLoader(
            self.data_set,
            batch_size=len(self.data_set),
            num_workers=3,
            drop_last=False, 
            shuffle=False
        )