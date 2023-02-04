import h5py
import torch
from pytorch_lightning import Trainer, seed_everything
from pytorch_lightning.loggers import TensorBoardLogger
from pytorch_lightning.callbacks import ModelCheckpoint, TQDMProgressBar
from pytorch_lightning.callbacks.early_stopping import EarlyStopping
import numpy as np
from modules.model import dab
from modules.wae import wae
from modules.data_loader import DataModule
from argparse import ArgumentParser
seed_everything(42, workers=True)
def train(hparams):
    train_loader = DataModule(data_name=hparams.dataset_name, batchsize=4096)
    train_loader.prepare_data()
    train_loader.setup()

    in_dim = train_loader.ann_data.n_vars
    dimentions = [in_dim, 512, 128, 32]

    # load bulk information
    # path = '/home/suyanchi/project/dab/data/' + hparams.dataset_name + '_bulk.csv'
    # bulk_info = np.loadtxt(path, delimiter=',')
    # build model
    autoEncoder = wae(dimentions, hparams.dataset_name, hparams.alpha)
    # autoEncoder = dab(dimentions, hparams.dataset_name, hparams.alpha)
    print(autoEncoder)

    checkpoint_callback = ModelCheckpoint(
        save_top_k=1,
        verbose=False,
        filename=hparams.dataset_name + '-{epoch:02d}-{val_loss:.2f}'
        # monitor='train_loss',
        # mode='min'
    )

    logger = TensorBoardLogger('log', name=hparams.dataset_name)
    early_stop_callback = EarlyStopping(monitor="val_loss", mode="min")
    trainer = Trainer(
        logger=logger,
        default_root_dir='./log',
        callbacks=[checkpoint_callback, early_stop_callback, TQDMProgressBar()],
        log_every_n_steps=50,
        max_epochs=hparams.max_epochs,
        accelerator='gpu',
        devices=[1],
        auto_select_gpus=True
        # auto_lr_find = True,
        # auto_scale_batch_size=True
    )

    # train
    trainer.fit(autoEncoder, train_loader)

    # predict 
    results = trainer.predict(autoEncoder, train_loader.test_dataloader())

    # save results cell *gene
    # re = (torch.exp(results[0][0])-1).numpy()
    re = (results[0][0]).numpy()
    f = h5py.File('/home/suyanchi/project/dab/results/1M/bis_4096.h5', 'w')
    f['data'] = re
    f.close

    



def main():
    parser = ArgumentParser()
    parser.add_argument('--gpus', default="1")
    parser.add_argument('--max_epochs', default=1000)
    parser.add_argument('--dataset_name', default="Spleen")
    parser.add_argument('--alpha', default=1)
    args = parser.parse_args()

    train(args)


if __name__ == "__main__":
    main()