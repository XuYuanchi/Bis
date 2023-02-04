import os
import anndata as ad
import numpy as np
import pandas as pd
from pandas import DataFrame, Series
import pytorch_lightning as pl
import scanpy as sc
from torch.utils.data import Dataset, DataLoader, random_split


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

    def prepare_data(self, filter_min_counts=False, size_factors=True, normalize_input=True, logtrans_input=True):
        # loading and proprecessing data
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
        )

    def val_dataloader(self):
        return DataLoader(
            self.dataset_val,
            batch_size=self.batch_size,
        )

    def test_dataloader(self):
        return DataLoader(
            self.dataset_test,
            batch_size=len(self.dataset_test)
        )