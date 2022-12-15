#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  7 15:18:03 2022

@author: omidard
"""
import cobra
import pandas as pd
from cobra.io import load_json_model
from glob import glob
from cobra.manipulation.delete import delete_model_genes, remove_genes
import os
from os.path import isfile, join
from Bio import Entrez, SeqIO
from os import listdir
from cobra.flux_analysis import gapfill         
import cobra
import pandas as pd
from glob import glob



re300_400=[]
re401_500=[]
re501_600=[]
re601_700=[]
re701_800=[]
re801_900=[]
re901_1000=[]
re1001_1100=[]
re1101_1200=[]
re1201_1300=[]
re1300_1400=[]
re1400_1500=[]
re1500_1600=[]
re1600_1700=[]
models2 =glob('%s/*.json'%'output_models_dir')
for ml in models2:
    model = load_json_model(ml)
    if len(model.reactions) > 300 and len(model.reactions) < 400:
        re300_400.append(model.id)
    if len(model.reactions) > 400 and len(model.reactions) < 500:
        re401_500.append(model.id)
    if len(model.reactions) > 500 and len(model.reactions) < 600:
        re501_600.append(model.id)
    if len(model.reactions) > 600 and len(model.reactions) < 700:
        re601_700.append(model.id)
    if len(model.reactions) > 700 and len(model.reactions) < 800:
        re701_800.append(model.id)
    if len(model.reactions) > 800 and len(model.reactions) < 900:
        re801_900.append(model.id)
    if len(model.reactions) > 900 and len(model.reactions) < 1000:
        re901_1000.append(model.id)
    if len(model.reactions) > 1000 and len(model.reactions) < 1100:
        re1001_1100.append(model.id)
    if len(model.reactions) > 1100 and len(model.reactions) < 1200:
        re1101_1200.append(model.id)
    if len(model.reactions) > 1200 and len(model.reactions) < 1300:
        re1201_1300.append(model.id)
    if len(model.reactions) > 1300 and len(model.reactions) < 1400:
        re1300_1400.append(model.id)
    if len(model.reactions) > 1400 and len(model.reactions) < 1500:
        re1400_1500.append(model.id)
    if len(model.reactions) > 1500 and len(model.reactions) < 1600:
        re1500_1600.append(model.id)
    if len(model.reactions) > 1600 and len(model.reactions) < 1700:
        re1600_1700.append(model.id)


#%%


model_id =[]
model_size =[]
models2 =glob('%s/*.json'%'output_models_dir')
for ml in models2:
    model = load_json_model(ml)
    model_id.append(model.id)
    model_size.append(len(model.reactions))
models3 = pd.DataFrame({'model_id':model_id,'model_size':model_size})





















