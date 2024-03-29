import numpy as np
import torch
from torch import nn
from torch import optim
import torch.nn.functional as F
import pytorch_lightning as pl
from .utils import js_divergence
from torchmetrics.functional import pearson_corrcoef

class MeanAct(nn.Module):
    def __init__(self):
        super(MeanAct, self).__init__()

    def forward(self, x):
        return torch.clamp(torch.exp(x), min=1e-5, max=1e6)

class DispAct(nn.Module):
    def __init__(self):
        super(DispAct, self).__init__()

    def forward(self, x):
        return torch.clamp(F.softplus(x), min=1e-4, max=1e4)

class ExpAct(nn.Module):
    def __init__(self):
        super(ExpAct, self).__init__()

    def forward(self, x):
        return torch.clamp(torch.exp(x), min=1e-8, max=1e8)

def get_zinb_loss(x,  mean, disp, pi,scale_factor=1.0, ridge_lambda=0.0):

    eps = 1e-10
    if isinstance(scale_factor,float):
        scale_factor=np.full((len(mean),),scale_factor)
    scale_factor = scale_factor[:, None]
    mean = mean * scale_factor

    t1 = torch.lgamma(disp+eps) + torch.lgamma(x+1.0) - torch.lgamma(x+disp+eps)
    t2 = (disp+x) * torch.log(1.0 + (mean/(disp+eps))) + (x * (torch.log(disp+eps) - torch.log(mean+eps)))
    nb_final = t1 + t2

    nb_case = nb_final - torch.log(1.0-pi+eps)
    zero_nb = torch.pow(disp/(disp+mean+eps), disp)
    zero_case = -torch.log(pi + ((1.0-pi)*zero_nb)+eps)
    result = torch.where(torch.le(x, 1e-8), zero_case, nb_case)

    if ridge_lambda > 0:
        ridge = ridge_lambda*torch.square(pi)
        result += ridge
    result = torch.mean(result)
    return result

class LAN(nn.Module):
    """
    Linear-Activation-Normalization.
    """

    def __init__(self, in_dim, out_dim, activation=nn.ReLU(True), norm=nn.Identity()):
        super(LAN, self).__init__()
        self.L = nn.Linear(in_dim, out_dim)
        self.B = nn.BatchNorm1d(out_dim)
        self.A = activation
        # self.N = norm

    def forward(self, x):
        z = self.L(x)
        z = self.B(z)
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

class dab(pl.LightningModule):
    def __init__(self, in_dimentions, data_name, alpha):
        super(dab, self).__init__()
        out_dimentions = list(reversed(in_dimentions))
        self.encoder = DENcoder(in_dimentions, nn.ReLU(True), nn.Identity())
        # self.decoder = DENcoder(out_dimentions, nn.ReLU(True), nn.Identity())
        self.decoder = DENcoder(out_dimentions[:-1], nn.ReLU(True), nn.Identity())
        self.dec_mean = nn.Sequential(nn.Linear(out_dimentions[-2], out_dimentions[-1]), MeanAct())
        self.dec_disp = nn.Sequential(nn.Linear(out_dimentions[-2], out_dimentions[-1]), DispAct())
        self.dec_pi = nn.Sequential(nn.Linear(out_dimentions[-2], out_dimentions[-1]), nn.Sigmoid())
        #input_dim = in_dimentions[0]
        #self.low_rank_layer = DENcoder([input_dim, input_dim, input_dim] ,nn.ReLU(False), nn.Identity())
        self.criterion = nn.MSELoss()
        bulk_vec = self.get_Bulk_data(data_name)
        self.register_buffer("bulk_vec", bulk_vec)
        self.alpha = alpha

    # get bulk data
    def get_Bulk_data(self, data_name):
        path = '/home/suyanchi/project/dab/data/' + data_name + '_bulk.csv'
        bulk_info = torch.from_numpy(np.loadtxt(path, delimiter=','))
        return torch.log1p(bulk_info).unsqueeze(0)

    def get_loss(self, x, x_hat, bulk_vec, alpha):
        mse_loss = nn.functional.mse_loss(x_hat, x)
        x_mean = torch.mean(x_hat, 0).unsqueeze(0)
        # r_loss = nn.functional.mse_loss(x_mean, bulk_vec)/pearson_corrcoef(x_mean, bulk_vec)
        r_loss = nn.functional.mse_loss(x_mean, bulk_vec)
        # r_loss = pearson_corrcoef(x_mean, bulk_vec)
        # CE = nn.CrossEntropyLoss()
        # re_loss = CE(x_mean, bulk_vec)
        re_loss = nn.functional.kl_div(F.softmax(x_mean, dim=-1).log(), F.softmax(bulk_vec, dim=-1))
        # re_loss = nn.functional.kl_div(x_mean, bulk_vec)
        # re_loss = sinkhorn_normalized(x_mean, bulk_vec, 0.01, 512, 100)
        # re_loss = js_divergence(x_mean, bulk_vec)
        loss = mse_loss + 0.4*re_loss + 0.01*r_loss
        return loss

    def get_loss1(self, x, mean, disp, pi, bulk_vec, alpha, sf):
        zinb_loss = get_zinb_loss(x, mean, disp, pi, sf)
        x_mean = torch.mean(mean, 0)
        # r_loss = nn.functional.mse_loss(x_mean, bulk_vec)/pearson_corrcoef(x_mean, bulk_vec)
        r_loss = nn.functional.mse_loss(x_mean, bulk_vec)
        # r_loss = pearson_corrcoef(x_mean, bulk_vec)
        # CE = nn.CrossEntropyLoss()
        # re_loss = CE(x_mean, bulk_vec)
        re_loss = nn.functional.kl_div(F.softmax(x_mean, dim=-1).log(), F.softmax(bulk_vec, dim=-1), reduction='sum')
        # re_loss = nn.functional.kl_div(x_mean, bulk_vec)
        # re_loss = sinkhorn_normalized(x_mean, bulk_vec, 0.01, 512, 100)
        # re_loss = js_divergence(x_mean, bulk_vec)
        loss = zinb_loss + 0*re_loss + 0.0*r_loss
        return loss



    def forward(self, x):
        z = self.encoder(x)
        h = self.decoder(z)
        mean = torch.exp(self.dec_mean(h))
        disp = torch.exp(self.dec_disp(h))
        pi = self.dec_pi(h)
        
        # x_hat = self.low_rank_layer(x_hat)
        return z, mean, disp, pi

    def training_step(self, batch, batch_idx):
        # training_step defines the train loop.
        # it is independent of forward
        x, sf = batch
        z, mean, disp, pi = self.forward(x)
        # add bulk constrains loss
        # loss = nn.functional.mse_loss(x_hat, x) + nn.functional.mse_loss(torch.mean(x_hat, 0), self.bulk_vec.float())
        # loss = self.get_loss(x, x_hat, self.bulk_vec.float(), self.alpha)
        loss = self.get_loss1(x, mean, disp, pi, self.bulk_vec.float(), self.alpha, sf)
        # Logging to TensorBoard by default
        self.log("train_loss", loss, prog_bar=True)
        return loss

    def test_step(self, batch, batch_idx):
        # training_step defines the train loop.
        # it is independent of forward
        x, sf = batch
        z, mean, disp, pi = self.forward(x)
        # add bulk constrains loss
        # loss = nn.functional.mse_loss(x_hat, x) + nn.functional.mse_loss(torch.mean(x_hat, 0), self.bulk_vec.float())
        # loss = self.get_loss(x, x_hat, self.bulk_vec.float(), self.alpha)
        loss = self.get_loss1(x, mean, disp, pi, self.bulk_vec.float(), self.alpha, sf)
        # Logging to TensorBoard by default
        self.log("test_loss", loss, prog_bar=True)

    def validation_step(self, batch, batch_idx):
        # training_step defines the train loop.
        # it is independent of forward
        x, sf = batch
        z, mean, disp, pi = self.forward(x)
        # add bulk constrains loss
        # loss = nn.functional.mse_loss(x_hat, x) + nn.functional.mse_loss(torch.mean(x_hat, 0), self.bulk_vec.float())
        # loss = self.get_loss(x, x_hat, self.bulk_vec.float(), self.alpha)
        loss = self.get_loss1(x, mean, disp, pi, self.bulk_vec.float(), self.alpha, sf)
        # Logging to TensorBoard by default
        self.log("val_loss", loss, prog_bar=True)
    
    def predict_step(self, batch, batch_idx, dataloader_idx=0):
        x, _ =batch
        return self(x)

    def configure_optimizers(self):
        optimizer = optim.Adam(self.parameters(), lr=1e-3, weight_decay=1e-5)
        return optimizer
