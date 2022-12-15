#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 09:44:46 2022

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
#%%
notes =glob('%s/*.txt'%'gap_inf_dir')
allinf = []
for tx in notes:
    with open(tx,'r') as txr:
        rx = set()
        for line in txr:
            line = line.strip()
            if line != "":
                cor= line.replace("', '","\n")
                cor2 = cor.replace("['","")
                cor3 =cor2.replace("']","")
                rx.add(cor3)
        allinf.append(rx)
coll =[]        
for i in range(len(allinf)):
    for z in allinf[i]:
        splitted = [x for x in z.split() if x != '']
        coll.append(splitted)
#%%

m11 = pd.DataFrame()
m6 = m1+m2+m3+m4+m5
m7 = [i for j, i in enumerate(m6) if i not in m6[:j]]


gapadd=[]
data_dir = 'ref_model_dir'
universal = cobra.io.read_sbml_model(join(data_dir, "marlbr2.xml"))
for reactions in universal.reactions:
    if reactions.id in m7:
        gapadd.append(reactions)

models =glob('%s/*.json'%'output_models_dir')
for ml in models:
    model = load_json_model(ml)
    mre =[]
    for re in model.reactions:
        mre.append(re.id)
    final = []
    for reactions in gapadd:
        if reactions.id not in mre:
            final.append(reactions)
    model.add_reactions(final)
    cobra.io.save_json_model(model, 'prep_models_dir/'+model.id)
            