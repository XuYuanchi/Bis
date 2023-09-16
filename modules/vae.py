import numpy as np
import torch
from torch import nn
from torch import optim
import pytorch_lightning as pl


class LAN(nn.Module):
    """
    Linear-Activation-Normalization.
    """

    def __init__(self, in_dim, out_dim, activation=nn.ReLU(True), norm=nn.Identity()):
        super(LAN, self).__init__()
        self.L = nn.Linear(in_dim, out_dim)
        self.A = activation
        self.N = norm

    def forward(self, x):
        z = self.L(x)
        z = self.A(z)
        z = self.N(z)
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

class dab_vae(pl.LightningModule):
    def __init__(self, in_dimentions, data_name, alpha):
        super(dab_vae, self).__init__()
        input_dim = in_dimentions[:-1]
        out_dimentions = list(reversed(in_dimentions))
        self.encoder = DENcoder(input_dim, nn.ReLU(True), nn.Identity())
        self.decoder = DENcoder(out_dimentions, nn.ReLU(True), nn.Identity())

        self.hidden2mu = nn.Linear(in_dimentions[-2], in_dimentions[-1])
        self.hidden2log_var = nn.Linear(in_dimentions[-2], in_dimentions[-1])
        self.log_scale = nn.Parameter(torch.Tensor([0.0]))
        #input_dim = in_dimentions[0]
        #self.low_rank_layer = DENcoder([input_dim, input_dim, input_dim] ,nn.ReLU(False), nn.Identity())
        self.criterion = nn.MSELoss()
        bulk_vec = self.get_Bulk_data(data_name)
        self.register_buffer("bulk_vec", bulk_vec)
        self.alpha = alpha
        self.kl_coeff = 1e-1

    # get bulk data
    def get_Bulk_data(self, data_name):
        path = '/home/suyanchi/project/dab/data/' + data_name + '_bulk.csv'
        bulk_info = torch.from_numpy(np.loadtxt(path, delimiter=','))
        return torch.log1p(bulk_info)

    def get_loss(self, x, x_hat, bulk_vec, alpha):
        x_mean = torch.mean(x_hat, 0)
        # r_loss = nn.functional.mse_loss(x_mean, bulk_vec)/pearson_corrcoef(x_mean, bulk_vec)
        r_loss = nn.functional.mse_loss(x_mean, bulk_vec)
        kl_loss = nn.functional.kl_div(x_mean, bulk_vec)
        loss = alpha*kl_loss + 0*r_loss
        return loss

    def gaussian_likelihood(self, x_hat, logscale, x):
        scale = torch.exp(logscale)
        mean = x_hat
        dist = torch.distributions.Normal(mean, scale)

        # measure prob of seeing image under p(x|z)
        log_pxz = dist.log_prob(x)

        return log_pxz.sum(dim=(1))

    def forward(self, x):
        z = self.encoder(x)
        mu = self.hidden2mu(z)
        log_var = self.hidden2log_var(z)

        std = torch.exp(log_var / 2)
        p = torch.distributions.Normal(torch.zeros_like(mu), torch.ones_like(std))
        q = torch.distributions.Normal(mu, std)
        z = q.rsample()

        #Push sample through decoder
        x_hat = self.decoder(z)

        return p, q, z, x_hat

    def training_step(self, batch, batch_idx):
        # training_step defines the train loop.
        # it is independent of forward
        x = batch
        p, q, z, x_hat = self.forward(x)
        # reconstruction loss
        recon_loss = nn.functional.mse_loss(x_hat, x, reduction="mean")

        #expectation under z of the kl divergence between q(z|x) and
        #a standard normal distribution of the same shape
        kl = torch.distributions.kl_divergence(q, p)
        kl = kl.mean()
        kl *= self.kl_coeff

        loss = kl + recon_loss
        # add bulk constrains loss
        # loss = nn.functional.mse_loss(x_hat, x) + nn.functional.mse_loss(torch.mean(x_hat, 0), self.bulk_vec.float())
        # loss = loss + self.get_loss(x, x_hat, self.bulk_vec.float(), self.alpha)
        # Logging to TensorBoard by default
        self.log("train_loss", loss, prog_bar=True)
        return loss

    def test_step(self, batch, batch_idx):
        # training_step defines the train loop.
        # it is independent of forward
        x = batch
        p, q, z, x_hat = self.forward(x)
        # reconstruction loss
        recon_loss = nn.functional.mse_loss(x_hat, x, reduction="mean")

        #expectation under z of the kl divergence between q(z|x) and
        #a standard normal distribution of the same shape
        kl = torch.distributions.kl_divergence(q, p)
        kl = kl.mean()
        kl *= self.kl_coeff

        loss = kl + recon_loss
        # add bulk constrains loss
        # loss = nn.functional.mse_loss(x_hat, x) + nn.functional.mse_loss(torch.mean(x_hat, 0), self.bulk_vec.float())
        # loss = loss + self.get_loss(x, x_hat, self.bulk_vec.float(), self.alpha)
        # Logging to TensorBoard by default
        self.log("test_loss", loss, prog_bar=True)
        return loss

    def validation_step(self, batch, batch_idx):
        # training_step defines the train loop.
        # it is independent of forward
        x = batch
        p, q, z, x_hat = self.forward(x)
        # reconstruction loss
        recon_loss = nn.functional.mse_loss(x_hat, x, reduction="mean")

        #expectation under z of the kl divergence between q(z|x) and
        #a standard normal distribution of the same shape
        kl = torch.distributions.kl_divergence(q, p)
        kl = kl.mean()
        kl *= self.kl_coeff

        loss = kl + recon_loss
        # add bulk constrains loss
        # loss = nn.functional.mse_loss(x_hat, x) + nn.functional.mse_loss(torch.mean(x_hat, 0), self.bulk_vec.float())
        # loss = loss + self.get_loss(x, x_hat, self.bulk_vec.float(), self.alpha)
        # Logging to TensorBoard by default
        self.log("val_loss", loss, prog_bar=True)
        return loss
    
    def predict_step(self, batch, batch_idx, dataloader_idx=0):
        return self(batch)

    def configure_optimizers(self):
        optimizer = optim.Adam(self.parameters(), lr=1e-3, weight_decay=1e-5)
        return optimizer
