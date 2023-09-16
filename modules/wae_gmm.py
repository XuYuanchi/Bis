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
import math
from sklearn.mixture import GaussianMixture


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
    def __init__(self, in_dimentions, data_name, alpha, n_centroids):
        super(wae, self).__init__()
        out_dimentions = list(reversed(in_dimentions))
        self.encoder = DENcoder(in_dimentions[:-1], nn.ReLU(True), nn.Identity())
        self.decoder = DENcoder(out_dimentions, nn.ReLU(True), nn.Identity())
        # self.drop = nn.Dropout(p=0.1)
        self.fc_mean = nn.Linear(in_dimentions[-2], in_dimentions[-1])
        self.fc_logvar = nn.Linear(in_dimentions[-2], in_dimentions[-1])
        # input_dim = in_dimentions[0]
        #self.low_rank_layer = DENcoder([input_dim, input_dim, input_dim] ,nn.ReLU(False), nn.Identity())
        # self.criterion = nn.MSELoss()

        self.n_centroids = n_centroids
        # init c_params
        self.pi = nn.Parameter(torch.ones(n_centroids)/n_centroids)  # pc
        self.mu_c = nn.Parameter(torch.zeros(in_dimentions[-1], n_centroids)) # mu
        self.var_c = nn.Parameter(torch.ones(in_dimentions[-1], n_centroids)) # sigma^2

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
        else:
            path = '/home/suyanchi/project/dab/data/' + data_name + '_bulk.csv'
            bulk_info = torch.from_numpy(np.loadtxt(path, delimiter=','))
        return torch.log1p(bulk_info).unsqueeze(0)

    def get_loss(self, x, x_hat, bulk_vec, alpha, gamma, c_params, z_params):

        # reconstruction loss
        # mse_loss = nn.functional.smooth_l1_loss(x_hat, x)
        mse_loss = nn.functional.mse_loss(x_hat, x)
        # kl loss
        mu, logvar = z_params
        # log q(z|x) or q entropy    
        qentropy = -0.5*torch.sum(1+logvar+math.log(2*math.pi), 1)
        # log q(c|x)
        logqcx = torch.sum(gamma*torch.log(gamma), 1)
        # log p(z|c)
        mu_c, var_c, pi = c_params; #print(mu_c.size(), var_c.size(), pi.size())
        var_c += 1e-8
        n_centroids = pi.size(1)
        mu_expand = mu.unsqueeze(2).expand(mu.size(0), mu.size(1), n_centroids)
        logvar_expand = logvar.unsqueeze(2).expand(logvar.size(0), logvar.size(1), n_centroids)
        logpzc = -0.5*torch.sum(gamma*torch.sum(math.log(2*math.pi) + torch.log(var_c) + torch.exp(logvar_expand)/var_c + (mu_expand-mu_c)**2/var_c, dim=1), dim=1)

        # log p(c)
        logpc = torch.sum(gamma*torch.log(pi), 1)
        # KL(q(z,c|x)||p(z,c))
        # - log p(z|c) - log p(c) + log q(z|x) + log q(c|x)
        kld = (-logpzc - logpc + qentropy + logqcx).mean()

        # p_z = torch.rand_like(q_z)
        # d_loss = mmd(q_z, p_z, "IMQ", 1)
        x_mean = torch.mean(x_hat, 0).unsqueeze(0)
        # r_loss = nn.functional.mse_loss(x_mean, bulk_vec)/pearson_corrcoef(x_mean, bulk_vec)
        # r_loss = nn.functional.mse_loss(x_mean, bulk_vec)
        # r_loss = pearson_corrcoef(x_mean, bulk_vec)
        # CE = nn.CrossEntropyLoss()
        # re_loss = CE(x_mean, bulk_vec)
        # re_loss = nn.functional.kl_div(F.softmax(x_mean, dim=-1).log(), F.softmax(bulk_vec, dim=-1))
        re_loss = nn.functional.kl_div(x_mean, bulk_vec)
        # re_loss = sinkhorn_normalized(x_mean, bulk_vec, 0.01, 512, 100)
        # re_loss = js_divergence(x_mean, bulk_vec)
        # re_loss = self.imq_kernel(x_mean, bulk_vec,x_hat.shape[1])
        loss = mse_loss + 0.00*kld + 0.3*re_loss
        return loss


    def get_gamma(self, z):
        """
        Inference c from z
        gamma is q(c|x)
        q(c|x) = p(c|z) = p(c)p(c|z)/p(z)
        """
        n_centroids = self.n_centroids

        N = z.size(0)
        z = z.unsqueeze(2).expand(z.size(0), z.size(1), n_centroids)
        pi = self.pi.repeat(N, 1) # NxK
#         pi = torch.clamp(self.pi.repeat(N,1), 1e-10, 1) # NxK
        mu_c = self.mu_c.repeat(N,1,1) # NxDxK
        var_c = self.var_c.repeat(N,1,1) + 1e-8 # NxDxK

        # p(c,z) = p(c)*p(z|c) as p_c_z
        p_c_z = torch.exp(torch.log(pi) - torch.sum(0.5*torch.log(2*math.pi*var_c) + (z-mu_c)**2/(2*var_c), dim=1)) + 1e-10
        gamma = p_c_z / torch.sum(p_c_z, dim=1, keepdim=True)

        return gamma, mu_c, var_c, pi

    def encodeBatch(self, dataloader, out='z'):
        output = []
        for x, sf in dataloader:
            x = x.view(x.size(0), -1).float()
            h = self.encoder(x)
            mean = self.fc_mean(h)
            logvar = self.fc_logvar(h)
            z = self.reparameterize(mean, logvar)

            if out == 'z':
                output.append(z.detach())
            elif out == 'x':
                recon_x = self.decoder(z)
                output.append(recon_x.data)
            elif out == 'logit':
                output.append(self.get_gamma(z)[0].data)

        output = torch.cat(output).numpy()

        return output
    def init_gmm_params(self, batch):
        gmm = GaussianMixture(n_components=self.n_centroids, covariance_type='diag')
        z = self.encodeBatch(batch)
        gmm.fit(z)
        self.mu_c.data.copy_(torch.from_numpy(gmm.means_.T.astype(np.float32)))
        self.var_c.data.copy_(torch.from_numpy(gmm.covariances_.T.astype(np.float32)))

    def reparameterize(self, mean, logvar):
        sd = torch.exp(0.5 * logvar)  # Standard deviation
        # We'll assume the posterior is a multivariate Gaussian
        eps = torch.randn_like(sd)
        z = eps.mul(sd).add(mean)
        return z

    def forward(self, x):
        h = self.encoder(x)
        mean = self.fc_mean(h)
        logvar = self.fc_logvar(h)
        z = self.reparameterize(mean, logvar)
        x_hat = self.decoder(z)
        # x_hat = self.low_rank_layer(x_hat)
        return x_hat, z, mean, logvar

    

    def training_step(self, batch, batch_idx):
        # training_step defines the train loop.
        # it is independent of forward
        x, sf = batch
        # x_drop = self.drop(x)
        x_hat, z, mean, logvar = self.forward(x)
        # add bulk constrains loss
        # loss = nn.functional.mse_loss(x_hat, x) + nn.functional.mse_loss(torch.mean(x_hat, 0), self.bulk_vec.float())
        # loss = self.get_loss(x, x_hat, self.bulk_vec.float(), self.alpha)
        gamma, mu_c, var_c, pi = self.get_gamma(z) #, self.n_centroids, c_params)
        loss = self.get_loss(x, x_hat, self.bulk_vec.float(), self.alpha, gamma, (mu_c, var_c, pi), (mean, logvar))
        # Logging to TensorBoard by default
        self.log("train_loss", loss, prog_bar=True)
        return loss

    def test_step(self, batch, batch_idx):
        # training_step defines the train loop.
        # it is independent of forward
        x, sf = batch
        x_hat, z, mean, logvar = self.forward(x)
        # add bulk constrains loss
        # loss = nn.functional.mse_loss(x_hat, x) + nn.functional.mse_loss(torch.mean(x_hat, 0), self.bulk_vec.float())
        # loss = self.get_loss(x, x_hat, self.bulk_vec.float(), self.alpha)
        gamma, mu_c, var_c, pi = self.get_gamma(z) #, self.n_centroids, c_params)
        loss = self.get_loss(x, x_hat, self.bulk_vec.float(), self.alpha, gamma, (mu_c, var_c, pi), (mean, logvar))
        # Logging to TensorBoard by default
        self.log("test_loss", loss, prog_bar=True)

    def validation_step(self, batch, batch_idx):
        # training_step defines the train loop.
        # it is independent of forward
        x, sf = batch
        x_hat, z, mean, logvar = self.forward(x)
        # add bulk constrains loss
        # loss = nn.functional.mse_loss(x_hat, x) + nn.functional.mse_loss(torch.mean(x_hat, 0), self.bulk_vec.float())
        # loss = self.get_loss(x, x_hat, self.bulk_vec.float(), self.alpha)
        gamma, mu_c, var_c, pi = self.get_gamma(z) #, self.n_centroids, c_params)
        loss = self.get_loss(x, x_hat, self.bulk_vec.float(), self.alpha, gamma, (mu_c, var_c, pi), (mean, logvar))
        # Logging to TensorBoard by default
        self.log("val_loss", loss, prog_bar=True)
    
    def predict_step(self, batch, batch_idx, dataloader_idx=0):
        x, _ =batch
        return self(x)
    def predict(self, data_loader):
        return self.encodeBatch(data_loader, out='x')

    def configure_optimizers(self):
        optimizer = optim.Adam(self.parameters(), lr=1e-3, weight_decay=1e-5)
        return optimizer
