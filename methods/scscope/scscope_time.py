import methods.scscope.scscope_ti as DeepImpute
import pandas as pd
import numpy as np
import scanpy as sc
import h5py
import time
def RUN_MAIN():

    # 1. Load gene expression matrix of simulated data
    # gene expression with simulated dropouts
    f = h5py.File('/home/suyanchi/project/dab/data/1M/1M.h5')
    # f = h5py.File('/home/suyanchi/project/dab/data/time/100000.h5')
    dat_sc = f['data']
    f.close
    t1 = time.time()
    dat_sc = pd.DataFrame(np.array(dat_sc)[:, 1:1000])

    # 2. Normalize gene expression based on scanpy (normalize each cell to have same library size)
    # matrix of cells x genes
    raw_library_size = dat_sc.sum(0)

    gene_expression = sc.AnnData(dat_sc.values)
    # normalize each cell to have same count number
    sc.pp.normalize_per_cell(gene_expression)
    # update datastructure to use normalized data
    gene_expression = gene_expression.X

    latent_dim = 50

    # 3. scScope learning
    if gene_expression.shape[0] >= 100000:
        DI_model = DeepImpute.train(gene_expression, latent_dim, T=2, batch_size=64, max_epoch=300, num_gpus=1)
    else:
        DI_model = DeepImpute.train(gene_expression, latent_dim, T=2, batch_size=64, max_epoch=300, num_gpus=1)

    # 4. latent representations and imputed expressions
    _, imputed_val, _ = DeepImpute.predict(gene_expression, DI_model)
    # print(imputed_val.shape)
    # print(raw_library_size.values.shape)
    # print(gene_expression.sum(1, keepdims=True).shape)
    # imputed_val_resume = imputed_val * raw_library_size / gene_expression.sum(1, keepdims=True)
    t2 = time.time()
    print(t2-t1)


if __name__ == '__main__':
    RUN_MAIN()