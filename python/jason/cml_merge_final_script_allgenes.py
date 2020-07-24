#!/usr/bin/env python
# coding: utf-8
# %%

# # Using scVI with the optimized parameters as in the examples in the scVI docs

# %%


import matplotlib
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns

import numpy as np
import numpy.random as random
import pandas as pd
import scanpy as sc
import louvain

use_cuda = True
from scvi.dataset.dataset import GeneExpressionDataset
from scvi.inference import UnsupervisedTrainer
from scvi.models import SCANVI, VAE
from scvi.dataset.csv import CsvDataset

from umap import UMAP


# %%


b2b = CsvDataset("batch2_b_qcmat.csv", save_path=pathway)
b26m = CsvDataset("batch2_6m_qcmat.csv", save_path=pathway)
b2r = CsvDataset("batch2_r_qcmat.csv", save_path=pathway)

b3b = CsvDataset("batch3_b_qcmat.csv", save_path=pathway)
b36m = CsvDataset("batch3_6m_qcmat.csv", save_path=pathway)
b3r = CsvDataset("batch3_r_qcmat.csv", save_path=pathway)

b4b = CsvDataset("batch4_b_qcmat.csv", save_path=pathway)
b412m = CsvDataset("batch4_12m_qcmat.csv", save_path=pathway)

b5b = CsvDataset("batch5_b_qcmat.csv", save_path=pathway)
b56m = CsvDataset("batch5_6m_qcmat.csv", save_path=pathway)
b512m = CsvDataset("batch5_12m_qcmat.csv", save_path=pathway)

b6b = CsvDataset("batch_6_b_mat.csv", save_path=pathway)
b6r = CsvDataset("batch_6_r_mat.csv", save_path=pathway)

b7b = CsvDataset("batch_7_b_mat.csv", save_path=pathway)
b7r = CsvDataset("batch_7_r_mat.csv", save_path=pathway)

all_dataset = GeneExpressionDataset()
all_dataset.populate_from_datasets([b2b, b26m, b2r, b3b, b36m, b3r, b4b, b412m, b5b, b56m, b512m, b6b, b6r, b7b, b7r])


# %%


vae = VAE(all_dataset.nb_genes, n_batch=all_dataset.n_batches, n_labels=all_dataset.n_labels,
          n_hidden=256, n_latent=12, n_layers=3, dropout_rate=0.1, reconstruction_loss='nb', dispersion='gene-cell')


# %%


trainer = UnsupervisedTrainer(vae, all_dataset, train_size=0.99)


# %%


trainer.train(n_epochs=200, lr=0.005)


# %%


full = trainer.create_posterior(trainer.model, all_dataset, indices=np.arange(len(all_dataset)))
latent, batch_indices, labels = full.sequential().get_latent()
batch_indices = batch_indices.ravel()


# %%


latent_u = UMAP(spread=2).fit_transform(latent)


# %%


#cm = LinearSegmentedColormap.from_list(
#        'my_cm', ['deepskyblue', 'hotpink'], N=2)
cm = plt.get_cmap('Paired') 
fig, ax = plt.subplots(figsize=(11, 8.5))
order = np.arange(latent.shape[0])
random.shuffle(order)
ax.scatter(latent_u[order, 0], latent_u[order, 1], 
           c=all_dataset.batch_indices.ravel()[order], 
           cmap=cm, edgecolors='none', s=1)    
plt.axis("off")
fig.set_tight_layout(True)
plt.savefig('umap_cml_final_11dec.pdf')


# %%


import torch

torch.save(trainer.model.state_dict(), 'cml_final_11dec.pt')
torch.save(trainer.model.state_dict(), 'cml_final_11dec.pkl')

np.savetxt("cml_final_latents_11dec.csv", latent, delimiter=",")
np.savetxt("cml_final_indices_11dec.csv", batch_indices, delimiter=",")

