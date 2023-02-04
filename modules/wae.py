import numpy as np
import torch
import h5py
import scanpy as sc
from anndata import AnnData
from torch.autograd import Variable
from torch import nn
from torch import optim
import torch.nn.functional as F
import pytorch_lightning as pl
from .utils import js_divergence
from torchmetrics.functional import pearson_corrcoef
import scipy.io as sio

def compute_kernel(x, y, kernel, imq_c):
    """
    gaussian kernel
    :param x: point x
    :param y: point y
    :return: kernel distance between x, y
    """
    x_size = x.shape[0]
    y_size = y.shape[0]
    dim = x.shape[1]

    tiled_x = x.view(x_size, 1, dim).repeat(1, y_size, 1)
    tiled_y = y.view(1, y_size, dim).repeat(x_size, 1, 1)
    euclidean_dist = torch.mean((tiled_x - tiled_y) ** 2, dim=2)
    if kernel == "RBF":
        computed_kernel = torch.exp(-euclidean_dist / dim * 1.0)
    elif kernel == "IMQ":
        computed_kernel = 1. / (euclidean_dist + imq_c)
    return computed_kernel

def mmd(x, y, kernel="RBF", imq_c=1):
    """
    mmd distance
    :param x: point x
    :param y: point y
    :return: mmd distance between x, y
    """
    x_kernel = compute_kernel(x, x, kernel, imq_c)
    y_kernel = compute_kernel(y, y, kernel, imq_c)
    xy_kernel = compute_kernel(x, y, kernel, imq_c)
    return torch.mean(x_kernel) + torch.mean(y_kernel) - 2 * torch.mean(xy_kernel)



class LAN(nn.Module):
    """
    Linear-Activation-Normalization.
    """

    def __init__(self, in_dim, out_dim, activation=nn.ReLU(True), norm=nn.Identity()):
        super(LAN, self).__init__()
        self.L = nn.Linear(in_dim, out_dim)
        #self.B = nn.BatchNorm1d(out_dim)
        self.A = activation
        # self.N = norm

    def forward(self, x):
        z = self.L(x)
        #z = self.B(z)
        z = self.A(z)
        # z = self.N(z)
        return z

class DENcoder(nn.Module):
    def __init__(self, dimentions, mid_activation, last_activation):
        super(DENcoder, self).__init__()
        layers = []
        in_dims = dimentions[:-1]

        for i, in_dim in enumerate(in_dims):
            # activation = last_activation if i + 1 == len(in_dims) else mid_activation
            out_dim = dimentions[i + 1]
            layers.append(
                LAN(in_dim, out_dim, activation=mid_activation)
            )

        self.model = nn.Sequential(*layers)

    def forward(self, x):
        return self.model(x)

class wae(pl.LightningModule):
    def __init__(self, in_dimentions, data_name, alpha):
        super(wae, self).__init__()
        out_dimentions = list(reversed(in_dimentions))
        self.encoder = DENcoder(in_dimentions, nn.ReLU(True), nn.Identity())
        self.decoder = DENcoder(out_dimentions, nn.ReLU(True), nn.Identity())
        bulk_vec = self.get_Bulk_data(data_name)
        self.register_buffer("bulk_vec", bulk_vec)
        self.alpha = alpha

    # get bulk data
    def get_Bulk_data(self, data_name):
        if data_name=='1M':
            f = h5py.File('/home/suyanchi/project/dab/data/1M/1M.h5')
            bulk = f['bulk']
            f.close
            bulk_adata = AnnData(np.array(bulk).T)
            sc.pp.normalize_total(bulk_adata, target_sum=1e6)
            bulk_info = torch.from_numpy(np.mean(bulk_adata.X, axis=1))
        elif data_name=='test':
            f = sio.loadmat('/home/suyanchi/project/dab/data/'+data_name+'/6/1.mat')
            #bulk_info = torch.from_numpy(np.array(f['data_true'], dtype=np.float32).T.mean())
            # bulk_info = torch.as_tensor(np.mean(np.array(f['data_true']).astype('float'), axis=1))
            bulk_adata = AnnData(np.array(f['data_true']).astype('float').T)
            sc.pp.normalize_total(bulk_adata, target_sum=1e6)
            bulk_info = torch.from_numpy(np.mean(bulk_adata.X, axis=0))
            # print(bulk_info.shape)
        elif data_name in ['liver', 'heart', 'marrow', 'lung']:
            f = sio.loadmat('/home/suyanchi/project/dab/data/downsample/'+data_name+'.mat')
            bulk_adata = AnnData(np.array(f['data_bulk']).astype('float').T)
            sc.pp.normalize_total(bulk_adata, target_sum=1e6)
            bulk_info = torch.from_numpy(np.mean(bulk_adata.X, axis=0))
        else:
            path = '/home/suyanchi/project/dab/data/' + data_name + '_bulk.csv'
            bulk_info = torch.from_numpy(np.loadtxt(path, delimiter=','))
        return torch.log1p(bulk_info).unsqueeze(0)

    def get_loss(self, x, x_hat, bulk_vec, alpha, q_z):

        # reconstruction loss
        # mse_loss = nn.functional.smooth_l1_loss(x_hat, x)
        mse_loss = nn.functional.mse_loss(x_hat, x)
        # d loss
        p_z = Variable(torch.rand_like(q_z)*1)
        mmd_loss = self.imq_kernel(q_z, p_z, h_dim=32)
        mmd_loss = mmd_loss / x_hat.shape[0]
        # p_z = torch.rand_like(q_z)
        # d_loss = mmd(q_z, p_z, "IMQ", 1)
        x_mean = torch.mean(x_hat, 0).unsqueeze(0)
        # r_loss = nn.functional.mse_loss(x_mean, bulk_vec)/pearson_corrcoef(x_mean, bulk_vec)
        r_loss = nn.functional.mse_loss(x_mean, bulk_vec)

        re_loss = nn.functional.kl_div(x_mean, bulk_vec)

        loss = mse_loss + 0.01*mmd_loss + 0.5*re_loss + 0*r_loss
        # clustering 0.01, 0.5
        return loss

    def imq_kernel(self, X: torch.Tensor, Y: torch.Tensor, h_dim: int):
        batch_size = X.size(0)

        norms_x = X.pow(2).sum(1, keepdim=True)  # batch_size x 1
        prods_x = torch.mm(X, X.t())  # batch_size x batch_size
        dists_x = norms_x + norms_x.t() - 2 * prods_x

        norms_y = Y.pow(2).sum(1, keepdim=True)  # batch_size x 1
        prods_y = torch.mm(Y, Y.t())  # batch_size x batch_size
        dists_y = norms_y + norms_y.t() - 2 * prods_y

        dot_prd = torch.mm(X, Y.t())
        dists_c = norms_x + norms_y.t() - 2 * dot_prd

        stats = 0
        for scale in [.1, .2, .5, 1., 2., 5., 10.]:
            C = 2 * h_dim * 1.0 * scale
            res1 = C / (C + dists_x)
            res1 += C / (C + dists_y)

            if torch.cuda.is_available():
                res1 = (1 - torch.eye(batch_size).to(self.device)) * res1
            else:
                res1 = (1 - torch.eye(batch_size)) * res1

            res1 = res1.sum() / (batch_size - 1)
            res2 = C / (C + dists_c)
            res2 = res2.sum() * 2. / (batch_size)
            stats += res1 - res2

        return stats

    def rbf_kernel(self, X: torch.Tensor, Y: torch.Tensor, h_dim: int):
        batch_size = X.size(0)

        norms_x = X.pow(2).sum(1, keepdim=True)  # batch_size x 1
        prods_x = torch.mm(X, X.t())  # batch_size x batch_size
        dists_x = norms_x + norms_x.t() - 2 * prods_x

        norms_y = Y.pow(2).sum(1, keepdim=True)  # batch_size x 1
        prods_y = torch.mm(Y, Y.t())  # batch_size x batch_size
        dists_y = norms_y + norms_y.t() - 2 * prods_y

        dot_prd = torch.mm(X, Y.t())
        dists_c = norms_x + norms_y.t() - 2 * dot_prd

        stats = 0
        for scale in [.1, .2, .5, 1., 2., 5., 10.]:
            C = 2 * h_dim * 1.0 / scale
            res1 = torch.exp(-C * dists_x)
            res1 += torch.exp(-C * dists_y)

            if torch.cuda.is_available():
                res1 = (1 - torch.eye(batch_size).to(self.device)) * res1
            else:
                res1 = (1 - torch.eye(batch_size)) * res1

            res1 = res1.sum() / (batch_size - 1)
            res2 = torch.exp(-C * dists_c)
            res2 = res2.sum() * 2. / batch_size
            stats += res1 - res2

        return stats



    def reparameterize(self, mean, logvar):
        sd = torch.exp(0.5 * logvar)  # Standard deviation
        # We'll assume the posterior is a multivariate Gaussian
        eps = torch.randn_like(sd)
        z = eps.mul(sd).add(mean)
        return z

    def forward(self, x):
        z = self.encoder(x)
        # mean = self.fc_mean(h)
        # logvar = self.fc_logvar(h)
        # z = self.reparameterize(mean, logvar)
        x_hat = self.decoder(z)
        # x_hat = self.low_rank_layer(x_hat)
        return x_hat, z

    def training_step(self, batch, batch_idx):
        # training_step defines the train loop.
        # it is independent of forward
        x, sf = batch
        # x_drop = self.drop(x)
        x_hat, q_z = self.forward(x)
        # add bulk constrains loss
        # loss = nn.functional.mse_loss(x_hat, x) + nn.functional.mse_loss(torch.mean(x_hat, 0), self.bulk_vec.float())
        # loss = self.get_loss(x, x_hat, self.bulk_vec.float(), self.alpha)
        loss = self.get_loss(x, x_hat, self.bulk_vec.float(), self.alpha, q_z)
        # Logging to TensorBoard by default
        self.log("train_loss", loss, prog_bar=True)
        return loss

    def test_step(self, batch, batch_idx):
        # training_step defines the train loop.
        # it is independent of forward
        x, sf = batch
        x_hat, q_z = self.forward(x)
        # add bulk constrains loss
        # loss = nn.functional.mse_loss(x_hat, x) + nn.functional.mse_loss(torch.mean(x_hat, 0), self.bulk_vec.float())
        # loss = self.get_loss(x, x_hat, self.bulk_vec.float(), self.alpha)
        loss = self.get_loss(x, x_hat, self.bulk_vec.float(), self.alpha, q_z)
        # Logging to TensorBoard by default
        self.log("test_loss", loss, prog_bar=True)

    def validation_step(self, batch, batch_idx):
        # training_step defines the train loop.
        # it is independent of forward
        x, sf = batch
        x_hat, q_z = self.forward(x)
        # add bulk constrains loss
        # loss = nn.functional.mse_loss(x_hat, x) + nn.functional.mse_loss(torch.mean(x_hat, 0), self.bulk_vec.float())
        # loss = self.get_loss(x, x_hat, self.bulk_vec.float(), self.alpha)
        loss = self.get_loss(x, x_hat, self.bulk_vec.float(), self.alpha, q_z)
        # Logging to TensorBoard by default
        self.log("val_loss", loss, prog_bar=True)
    
    def predict_step(self, batch, batch_idx, dataloader_idx=0):
        x, _ =batch
        return self(x)

    def configure_optimizers(self):
        optimizer = optim.Adam(self.parameters(), lr=1e-3, weight_decay=1e-5)
        return optimizer
