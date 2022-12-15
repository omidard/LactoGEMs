#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 13:03:32 2022

@author: omidard


collect biochemical information for generated models


"""
#%%imports
import pandas as pd
from cobra.io import load_json_model
from glob import glob
from cobra.manipulation.delete import delete_model_genes, remove_genes
import os
from os.path import isfile, join
from Bio import Entrez, SeqIO
from os import listdir
from cobra.flux_analysis import gapfill
from interruptingcow import timeout
import signal
from Bio import SeqIO
#%% collect model id for gapfilled models
models =glob('%s/*.json.json'%'gapfilled_models_dir')

modelid = []
for ml in models:
    model = load_json_model(ml)
    model.solver = "glpk"
    mi = model.id
    modelid.append(mi)
#%%collect metadata from related genomes
directory1 =r'/Users/omidard/Desktop/backupgenome'
files = glob('%s/*.gbk'%directory1)
files2 = sorted(files)
strain_info = []
for file in files2:
        handle = open(file)
        record = list(SeqIO.parse(handle, "genbank"))
        for i in record:
            for f in i.features:
                if f.type=='source':
                    info = {}
                    info['file'] = file.replace('/Users/omidard/Desktop/backupgenome/','')
                    info['id'] = file.replace('/Users/omidard/Desktop/backupgenome/','')
                    info['genome_size'] = len(i.seq)/10000
                    for q in f.qualifiers.keys():
                        info[q] = '|'.join(f.qualifiers[q])
                        strain_info.append(info)
sinf = pd.DataFrame(strain_info)
sinf2 = sinf.drop_duplicates(subset='id', keep='first', inplace=False, ignore_index=False)
sinf3 = sinf2.sort_values(by=['id'])
print(sinf3)

#%%sum up information for genomes which has a gapfilled model
modgen =[]
for i in modelid:
    modtogen = i.replace('.json','.gbk')
    modgen.append(modtogen)

modelinf =[]
strain=[]
organism=[]
for i in sinf3.index:
    x = sinf3.id[i]
    if x in modgen:
        modelinf.append(x)
        strain.append(sinf3.strain[i])
        organism.append(sinf3.organism[i])

modelsid=[]
for i in modelinf:
    json = i.replace('.gbk','.json.json')
    modelsid.append(json)


#FOLLOWING DATAFRAME CONTAINS INFORMATION ABOUT ORGANISM AND STRAIN OF EACH GAPFILLED MODELS 
models_information = pd.DataFrame({'models_id':modelsid,'organism':organism,'strain':strain})
#%%run FBA for each models and collect data upon secreted metabolites



        
        
        
    













