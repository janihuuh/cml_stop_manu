
import os
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



### Train the PR1 model with 10 CV, save the model

_,_,_,_ = tcrgp.loso('/Users/janihuuh/Dropbox/cml_stop/results/tcrgp/input/pr1_tcrgp.csv',
                    'human',
                    'pr1',
                    pc_blo, cv=10,
                    cdr_types=[[],['cdr3']],
                    m_iters=2000,lr=0.005,nZ=0,mbs=0,
                    check_v="deduce",
                    va=None,vb='vb',cdr3a=None,cdr3b='cdr3b',epis='epitope',subs='subject')
plt.savefig('/Users/janihuuh/Dropbox/cml_stop/results/tcrgp/input/pr1_tcrgp.pdf')

 auc, params = tcrgp.train_classifier('/Users/janihuuh/Dropbox/cml_stop/results/tcrgp/input/pr1_tcrgp.csv','human',"pr1",pc_blo,
                                     cdr_types=[[],['cdr3']],m_iters=5000,lr=0.005,nZ=0,mbs=0,
                                     va=None,vb='vb',cdr3a=None,cdr3b='cdr3b',epis='epitope')

 with open('/Users/janihuuh/Dropbox/cml_stop/results/tcrgp/input/pr1_tcrgp','wb') as f:
    pickle.dump(params,f)

 _,_,_,_ = tcrgp.loso('/Users/janihuuh/Dropbox/tcrgp_run/src/TCRGP/training_data/paper/vdj_human_ATDALMTGY.csv','human','ATDALMTGY',pc_blo,cdr_types=[[],['cdr3']], m_iters=500,lr=0.005,nZ=0,mbs=0,va='va',vb='vb',cdr3a=None,cdr3b='cdr3b',epis='epitope',subs='subject')



##### Open the previously made models
with open('/Users/janihuuh/Dropbox/cml_stop/results/tcrgp/input/pr1_tcrgp','rb') as f: pr1_cdr3b = pickle.load(f)
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


with open('models/tcra/model_dash_BMLF_cdr3a', 'rb') as f: BMLF_cdr3a = pickle.load(f)
with open('models/tcra/model_dash_BMLF_tra', 'rb') as f: BMLF_tra = pickle.load(f)
with open('models/tcra/model_dash_M1_cdr3a', 'rb') as f: M1_cdr3a = pickle.load(f)
with open('models/tcra/model_dash_M1_tra', 'rb') as f: M1_tra = pickle.load(f)
with open('models/tcra/model_dash_pp65_cdr3a', 'rb') as f: pp65_cdr3a = pickle.load(f)
with open('models/tcra/model_dash_pp65_tra', 'rb') as f: pp65_tra = pickle.load(f)

tra_models = [BMLF_cdr3a, BMLF_tra, M1_cdr3a, M1_tra, pp65_cdr3a, pp65_tra]
tra_model_names = ["BMLF_cdr3a", "BMLF_tra", "M1_cdr3a",  "M1_tra", "pp65_cdr3a", "pp65_tra"]

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


models = viral_models + melanoma_models
model_names = viral_model_names + melanoma_model_names





### CML predictions emerson batch 1/5
models = [pr1_cdr3b] + viral_models + melanoma_models
model_names = ["pr1_cdr3b"] + viral_model_names + melanoma_model_names
file="/Users/"+me+"/Dropbox/cml_stop/results/tcrgp/healthy_input1.txt"
filename="healthy1"
output_dir="/Users/"+me+"/Dropbox/cml_stop/results/tcrgp/raw/healthy/"

## Make predictions with each of the models
for model, model_name in zip(models, model_names):
    preds = tcrgp.predict(file, model, cdr3b="cdr3aa", vb="v", delimiter='\t')
    with open(output_dir + filename + "_" + model_name + ".csv", 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        writer.writerow(['prediction'])
        for p in preds:
            writer.writerow(['{:.4f}'.format(p[0])])

### CML predictions emerson  batch 2/5
filename="healthy2"
output_dir="/Users/"+me+"/Dropbox/cml_stop/results/tcrgp/raw/healthy/"

for model, model_name in zip(models, model_names):
    preds = tcrgp.predict(file, model, cdr3b="cdr3aa", vb="v", delimiter='\t')
    with open(output_dir + filename + "_" + model_name + ".csv", 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        writer.writerow(['prediction'])
        for p in preds:
            writer.writerow(['{:.4f}'.format(p[0])])

### CML predictions emerson  batch 3/5
file="/Users/"+me+"/Dropbox/cml_stop/results/tcrgp/healthy_input3.txt"
filename="healthy3"

for model, model_name in zip(models, model_names):
    preds = tcrgp.predict(file, model, cdr3b="cdr3aa", vb="v", delimiter='\t')
    with open(output_dir + filename + "_" + model_name + ".csv", 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        writer.writerow(['prediction'])
        for p in preds:
            writer.writerow(['{:.4f}'.format(p[0])])

### CML predictions emerson  batch 4/5
file="/Users/"+me+"/Dropbox/cml_stop/results/tcrgp/healthy_input4.txt"
filename="healthy4"

for model, model_name in zip(models, model_names):
    preds = tcrgp.predict(file, model, cdr3b="cdr3aa", vb="v", delimiter='\t')
    with open(output_dir + filename + "_" + model_name + ".csv", 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        writer.writerow(['prediction'])
        for p in preds:
            writer.writerow(['{:.4f}'.format(p[0])])

### CML predictions emerson  batch 5/5
file="/Users/"+me+"/Dropbox/cml_stop/results/tcrgp/healthy_input5.txt"
filename="healthy5"

for model, model_name in zip(models, model_names):
    preds = tcrgp.predict(file, model, cdr3b="cdr3aa", vb="v", delimiter='\t')
    with open(output_dir + filename + "_" + model_name + ".csv", 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        writer.writerow(['prediction'])
        for p in preds:
            writer.writerow(['{:.4f}'.format(p[0])])




####### CML-stop predictions
models = [pr1_cdr3b] + viral_models + melanoma_models
model_names = ["pr1_cdr3b"] + viral_model_names + melanoma_model_names

file="/Users/"+me+"/Dropbox/cml_stop/results/tcrgp/cml_input.txt"
filename="cml"
output_dir="/Users/"+me+"/Dropbox/cml_stop/results/tcrgp/raw/"

## Make predictions with each of the models
for model, model_name in zip(models, model_names):
  preds = tcrgp.predict(file, model, cdr3b="cdr3aa", vb="v", delimiter='\t')

## Write results
with open(output_dir + "tcrb/" + filename + "_" + model_name + ".csv", 'w', newline='') as csvfile:
  writer = csv.writer(csvfile, delimiter=',')
writer.writerow(['prediction'])
for p in preds:
  writer.writerow(['{:.4f}'.format(p[0])])
