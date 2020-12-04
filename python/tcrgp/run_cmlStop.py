import os
me="janihuuh"
os.chdir("/Users/janihuuh/Dropbox/MelanoMAP/src/jani/python/tcrgp/")

import tcrgp
import pickle
import ast
import csv
import numpy as np
from matplotlib import pyplot as plt
plt.style.use('fivethirtyeight')

import pandas as pd
from glob import glob

subsmat = tcrgp.subsmatFromAA2('HENS920102')
pc_blo = tcrgp.get_pcs(subsmat,d=21)

cdrs = tcrgp.create_cdr_dict(alignment='imgt',species=['human'])

with open('models/common_viral/model_GILGFVFTL_cdr3ab','rb') as f: GILGFVFTL_cdr3ab = pickle.load(f)

with open('models/common_viral/model_GILGFVFTL_cdr3b','rb') as f: GILGFVFTL_cdr3b = pickle.load(f)
with open('models/common_viral/model_GLCTLVAML_cdr3ab','rb') as f: GLCTLVAML_cdr3ab = pickle.load(f)
with open('models/common_viral/model_GLCTLVAML_cdr3b','rb') as f: GLCTLVAML_cdr3b = pickle.load(f)
with open('models/common_viral/model_IPSINVHHY_cdr3b','rb') as f: IPSINVHHY_cdr3b = pickle.load(f)
with open('models/common_viral/model_NLVPMVATV_cdr3ab','rb') as f: NLVPMVATV_cdr3ab = pickle.load(f)
with open('models/common_viral/model_NLVPMVATV_cdr3b','rb') as f: NLVPMVATV_cdr3b = pickle.load(f)
with open('models/common_viral/model_PKYVKQNTLKLAT_cdr3b','rb') as f: PKYVKQNTLKLAT_cdr3b = pickle.load(f)
with open('models/common_viral/model_RAKFKQLL_cdr3b','rb') as f: RAKFKQLL_cdr3b = pickle.load(f)
with open('models/common_viral/model_RPRGEVRFL_cdr3b','rb') as f: RPRGEVRFL_cdr3b = pickle.load(f)
with open('models/common_viral/model_TPRVTGGGAM_cdr3b','rb') as f: TPRVTGGGAM_cdr3b = pickle.load(f)
with open('models/common_viral/model_YVLDHLIVV_cdr3b','rb') as f: YVLDHLIVV_cdr3b = pickle.load(f)

with open('models/melanoma/model_AMFWSVPTV_cdr3b','rb')       as f: AMFWSVPTV_cdr3b = pickle.load(f)
with open('models/melanoma/model_AMFWSVPTV_trb','rb')         as f: AMFWSVPTV_trb = pickle.load(f)
with open('models/melanoma/model_ELAGIGILTV_cdr3a_comb','rb') as f: ELAGIGILTV_cdr3a_comb = pickle.load(f)
with open('models/melanoma/model_ELAGIGILTV_cdr3ab','rb')     as f: ELAGIGILTV_cdr3ab = pickle.load(f)
with open('models/melanoma/model_ELAGIGILTV_cdr3b_comb','rb') as f: ELAGIGILTV_cdr3b_comb = pickle.load(f)
with open('models/melanoma/model_ELAGIGILTV_tra_comb','rb')   as f: ELAGIGILTV_tra_comb = pickle.load(f)
with open('models/melanoma/model_ELAGIGILTV_trab','rb')       as f: ELAGIGILTV_trab = pickle.load(f)
with open('models/melanoma/model_ELAGIGILTV_trb_comb','rb')   as f: ELAGIGILTV_trb_comb = pickle.load(f)
with open('models/melanoma/model_FLYNLLTRV_cdr3b','rb')       as f: FLYNLLTRV_cdr3b = pickle.load(f)
with open('models/melanoma/model_FLYNLLTRV_trb','rb')         as f: FLYNLLTRV_trb = pickle.load(f)
with open('models/melanoma/model_mart1_cdr3b','rb')           as f: mart1_cdr3b = pickle.load(f)
with open('models/melanoma/model_mart1_trb','rb')             as f: mart1_trb = pickle.load(f)
with open('models/melanoma/model_melana_cdr3a','rb')          as f: melana_cdr3a = pickle.load(f)
with open('models/melanoma/model_melana_cdr3b','rb')          as f: melana_cdr3b = pickle.load(f)
with open('models/melanoma/model_meloe1_cdr3a','rb')          as f: meloe1_cdr3a = pickle.load(f)
with open('models/melanoma/model_meloe1_cdr3b','rb')          as f: meloe1_cdr3b = pickle.load(f)


## The models
viral_models = [GILGFVFTL_cdr3b, GLCTLVAML_cdr3b, IPSINVHHY_cdr3b,
                NLVPMVATV_cdr3b, PKYVKQNTLKLAT_cdr3b, RAKFKQLL_cdr3b,
               RPRGEVRFL_cdr3b, TPRVTGGGAM_cdr3b, YVLDHLIVV_cdr3b]
viral_model_names = ["GILGFVFTL_cdr3b", "GLCTLVAML_cdr3b",     "IPSINVHHY_cdr3b",
                     "NLVPMVATV_cdr3b", "PKYVKQNTLKLAT_cdr3b", "RAKFKQLL_cdr3b",
                     "RPRGEVRFL_cdr3b", "TPRVTGGGAM_cdr3b",    "YVLDHLIVV_cdr3b"]

melanoma_models = [ELAGIGILTV_cdr3b_comb, FLYNLLTRV_cdr3b, mart1_cdr3b,
          melana_cdr3b, meloe1_cdr3b, GILGFVFTL_cdr3b, GLCTLVAML_cdr3b, IPSINVHHY_cdr3b,
          NLVPMVATV_cdr3b, PKYVKQNTLKLAT_cdr3b, RAKFKQLL_cdr3b,
          RPRGEVRFL_cdr3b, TPRVTGGGAM_cdr3b, YVLDHLIVV_cdr3b]
melanoma_model_names = ["ELAGIGILTV_cdr3b_comb", "FLYNLLTRV_cdr3b", "mart1_cdr3b", "melana_cdr3b", "meloe1_cdr3b",
                       "GILGFVFTL_cdr3b", "GLCTLVAML_cdr3b",     "IPSINVHHY_cdr3b",
                     "NLVPMVATV_cdr3b", "PKYVKQNTLKLAT_cdr3b", "RAKFKQLL_cdr3b",
                     "RPRGEVRFL_cdr3b", "TPRVTGGGAM_cdr3b",    "YVLDHLIVV_cdr3b"]

models = viral_models + melanoma_models
model_names = viral_model_names + melanoma_model_names

# ####### CML-stop predictions
file="/Users/"+me+"/Dropbox/cml_stop/results/tcrgp/cml_input.txt"
filename="cml"
output_dir="/Users/"+me+"/Dropbox/cml_stop/results/tcrgp/raw/"

## Make predictions with each of the models
for model, model_name in zip(models, model_names):
    preds = tcrgp.predict(file,
                          model,
                          cdr3b="cdr3aa", vb="v", delimiter='\t')

    ## Write results
    with open(output_dir + filename + "_" + model_name + ".csv", 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        writer.writerow(['prediction'])
        for p in preds:
            writer.writerow(['{:.4f}'.format(p[0])])
