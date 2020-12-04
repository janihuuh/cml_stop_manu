## Init
import os
import numpy as np
import numpy.random as random
import pandas as pd

from scvi.dataset.dataset import GeneExpressionDataset
from scvi.dataset.csv import CsvDataset
from scvi.inference import UnsupervisedTrainer
from scvi.models import SCANVI, VAE
from scvi.inference.autotune import auto_tune_scvi_model

from umap import UMAP

import torch
import scanpy as sc
import louvain

import logging
import pickle
from hyperopt import hp


# %matplotlib inline

use_cuda = True
n_epochs_all = None
save_path = ''
show_plot = True
os.chdir("/scratch/cs/csb/projects/cml_stop/")


## Read in the files

## CML stop
Batch_2_baseline = CsvDataset(filename='results/scvi/input_files/Batch_2_baseline.csv', save_path='', sep=',', new_n_genes=False)
Batch_2_6m = CsvDataset(filename='results/scvi/input_files/Batch_2_6m.csv', save_path='', sep=',', new_n_genes=False)
Batch_2_relapse = CsvDataset(filename='results/scvi/input_files/Batch_2_relapse.csv', save_path='', sep=',', new_n_genes=False)
Batch_3_baseline = CsvDataset(filename='results/scvi/input_files/Batch_3_baseline.csv', save_path='', sep=',', new_n_genes=False)
Batch_3_6m = CsvDataset(filename='results/scvi/input_files/Batch_3_6m.csv', save_path='', sep=',', new_n_genes=False)
Batch_3_relapse = CsvDataset(filename='results/scvi/input_files/Batch_3_relapse.csv', save_path='', sep=',', new_n_genes=False)
Batch_4_baseline = CsvDataset(filename='results/scvi/input_files/Batch_4_baseline.csv', save_path='', sep=',', new_n_genes=False)
Batch_4_12m = CsvDataset(filename='results/scvi/input_files/Batch_4_12m.csv', save_path='', sep=',', new_n_genes=False)
Batch_5_baseline = CsvDataset(filename='results/scvi/input_files/Batch_5_baseline.csv', save_path='', sep=',', new_n_genes=False)
Batch_5_6m = CsvDataset(filename='results/scvi/input_files/Batch_5_6m.csv', save_path='', sep=',', new_n_genes=False)
Batch_5_12m = CsvDataset(filename='results/scvi/input_files/Batch_5_12m.csv', save_path='', sep=',', new_n_genes=False)
Batch_6_baseline = CsvDataset(filename='results/scvi/input_files/Batch_6_baseline.csv', save_path='', sep=',', new_n_genes=False)
Batch_6_relapse = CsvDataset(filename='results/scvi/input_files/Batch_6_relapse.csv', save_path='', sep=',', new_n_genes=False)
Batch_7_baseline = CsvDataset(filename='results/scvi/input_files/Batch_7_baseline.csv', save_path='', sep=',', new_n_genes=False)
Batch_7_relapse = CsvDataset(filename='results/scvi/input_files/Batch_7_relapse.csv', save_path='', sep=',', new_n_genes=False)

## Wu Nature 2020 (5')
LB6 = CsvDataset(filename='results/scvi/input_files/LB6.csv', save_path='', sep=',', new_n_genes=False)
RB1 = CsvDataset(filename='results/scvi/input_files/RB1.csv', save_path='', sep=',', new_n_genes=False)
RB2 = CsvDataset(filename='results/scvi/input_files/RB2.csv', save_path='', sep=',', new_n_genes=False)
RB3 = CsvDataset(filename='results/scvi/input_files/RB3.csv', save_path='', sep=',', new_n_genes=False)

## Reider Nat Commun 2020 (3')
CLL6_d0 = CsvDataset(filename='results/scvi/input_files/CLL6_d0.csv', save_path='', sep=',', new_n_genes=False)
CLL8_d0 = CsvDataset(filename='results/scvi/input_files/CLL8_d0.csv', save_path='', sep=',', new_n_genes=False)
CLL1_d0 = CsvDataset(filename='results/scvi/input_files/CLL1_d0.csv', save_path='', sep=',', new_n_genes=False)
CLL5_d0 = CsvDataset(filename='results/scvi/input_files/CLL5_d0.csv', save_path='', sep=',', new_n_genes=False)

## 10X healthy (5')
healthy_10x = CsvDataset(filename='results/scvi/input_files/healthy_10x.csv', save_path='', sep=',', new_n_genes=False)



all_dataset = GeneExpressionDataset()
all_dataset.populate_from_per_batch_list(Xs =  [Batch_2_baseline.X,
                                                Batch_2_6m.X,
                                                Batch_2_relapse.X,

                                                Batch_3_baseline.X,
                                                Batch_3_6m.X,
                                                Batch_3_relapse.X,

                                                Batch_4_baseline.X,
                                                Batch_4_12m.X,

                                                Batch_5_baseline.X,
                                                Batch_5_6m.X,
                                                Batch_5_12m.X,

                                                Batch_6_baseline.X,
                                                Batch_6_relapse.X,

                                                Batch_7_baseline.X,
                                                Batch_7_relapse.X,

                                                LB6.X,
                                                RB1.X,
                                                RB2.X,
                                                RB3.X,

                                                CLL6_d0.X,
                                                CLL8_d0.X,
                                                CLL1_d0.X,
                                                CLL5_d0.X,

                                                healthy_10x.X])



## Train, save and fin
vae      = VAE(all_dataset.nb_genes,
               n_batch=all_dataset.n_batches,
               n_labels=all_dataset.n_labels,
               n_hidden=128,
               n_latent=30,
               n_layers=2,
               dispersion='gene')
trainer  = UnsupervisedTrainer(vae, all_dataset, train_size=1.0)
trainer.train(n_epochs=50)
torch.save(trainer.model.state_dict(), 'results/scvi/output/public_oneshot.pkl')


## Sample posterior to get latent representation and save those embeddings
full = trainer.create_posterior(trainer.model, all_dataset, indices=np.arange(len(all_dataset)))

latent, batch_indices, labels = full.sequential().get_latent()
batch_indices = batch_indices.ravel()

np.savetxt("results/scvi/output/public_oneshot_latent.csv", latent, delimiter=",")
np.savetxt("results/scvi/output/public_oneshot_indices.csv", batch_indices, delimiter=",")
