#!/usr/bin/env python
# coding: utf-8
# %%

# # CML Autotune for all samples
# Similar as the example in scVI docs but with a wider range of parameters for thorough search.

# %%


import sys

import logging
import os
import pickle

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
from hyperopt import hp

from scvi.inference import UnsupervisedTrainer
from scvi.inference.autotune import auto_tune_scvi_model
from scvi.models import VAE
from scvi.dataset.csv import CsvDataset
from scvi.dataset.dataset import GeneExpressionDataset


# %%


b2b = CsvDataset("batch2_b_qcmat.csv", save_path=pathname)
b26m = CsvDataset("batch2_6m_qcmat.csv", save_path=pathname)
b2r = CsvDataset("batch2_r_qcmat.csv", save_path=pathname)

b3b = CsvDataset("batch3_b_qcmat.csv", save_path=pathname)
b36m = CsvDataset("batch3_6m_qcmat.csv", save_path=pathname)
b3r = CsvDataset("batch3_r_qcmat.csv", save_path=pathname)

b4b = CsvDataset("batch4_b_qcmat.csv", save_path=pathname)
b412m = CsvDataset("batch4_12m_qcmat.csv", save_path=pathname)

b5b = CsvDataset("batch5_b_qcmat.csv", save_path=pathname)
b56m = CsvDataset("batch5_6m_qcmat.csv", save_path=pathname)
b512m = CsvDataset("batch5_12m_qcmat.csv", save_path=pathname)

b6b = CsvDataset("batch_6_b_mat.csv", save_path=pathname)
b6r = CsvDataset("batch_6_r_mat.csv", save_path=pathname)

b7b = CsvDataset("batch_7_b_mat.csv", save_path=pathname)
b7r = CsvDataset("batch_7_r_mat.csv", save_path=pathname)

all_dataset = GeneExpressionDataset()
all_dataset.populate_from_datasets([b2b, b26m, b2r, b3b, b36m, b3r, b4b, b412m, b5b, b56m, b512m, b6b, b6r, b7b, b7r])


# %%

# %%


n_epochs = 28
max_evals = 100
reserve_timeout = 180
fmin_timeout = 300

logger = logging.getLogger("scvi.inference.autotune")
logger.setLevel(logging.WARNING)

space = {
    "model_tunable_kwargs": {"n_latent": hp.choice("n_latent", np.arange(1, 31, 1)), 
                             "n_hidden": hp.choice("n_hidden", [32, 64, 128, 256]),
#                             "n_batch": hp.choice("n_batch", [all_dataset.n_batches]),
#                             "n_labels": hp.choice("n_labels", [all_dataset.n_labels]),
                             "n_layers": hp.choice("n_layer", [1, 2, 3, 4, 5]),
                             "dropout_rate": hp.choice("dropout_rate", [0.1, 0.3, 0.5, 0.7]),
                             "reconstruction_loss": hp.choice("reconstruction_loss", ['nb', 'zinb']),
                             "dispersion": hp.choice("dispersion", ['gene', 'gene-batch', 'gene-cell'])},
    "train_func_tunable_kwargs": {"lr": hp.choice("lr", [0.01, 0.005, 0.001, 0.0005, 0.0001])},
}


# %%


best_trainer, trials = auto_tune_scvi_model(
    gene_dataset=all_dataset,
    parallel=False,
    exp_key="cml_all",
    space = space,
    train_func_specific_kwargs={"n_epochs": n_epochs},
    max_evals=max_evals,
    save_path = savepath,
    reserve_timeout=reserve_timeout,
    fmin_timeout=fmin_timeout,
    use_batches = True
)

