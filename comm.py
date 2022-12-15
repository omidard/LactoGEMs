#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 18:12:41 2022
micom
@author: omidard
"""
from micom import load_pickle
import pandas as pd
from micom.workflows import build
from micom.workflows import grow
from micom.workflows import tradeoff

data=pd.read_csv('/home/omidard/gems/gems2/micom.csv')
model_db = '/home/omidard/gems/gems2'
manifest = build(data,out_folder="/home/omidard/models",model_db=None,cutoff=0.0001,threads=1)
#print(manifest)
com = load_pickle("/home/omidard/models/sample_1.pickle")

ex_ids = [r.id for r in com.exchanges]
print(ex_ids)
#print(com)
#print(len(com.reactions))
medium= pd.read_csv('/home/omidard/gems/medium.csv')
medium.drop(['reaction.1'],axis=1,inplace=True)
#medium['reaction'] = [i+'_m' for i in medium['reaction']]
#medium['metabolite'] = [i+'_m' for i in medium['metabolite']]
medium.set_index('reaction',inplace=True)
medium['reaction']=[i for i in medium.index]
medium = medium[['reaction', 'flux', 'metabolite']]
print(medium)
#manifest= pd.read_csv('/home/omidard/models/manifest.csv')
#manifest['file']=['/home/omidard/models/'+manifest['file'] for i in manifest['file']]
res = grow(manifest, model_folder="/home/omidard/models", medium=medium, tradeoff=0.5, threads=2)
tradeoff_rates = tradeoff(manifest, model_folder="/home/omidard/models", medium=medium, threads=2)
tradeoff_rates.head()
tradeoff_rates.groupby("tradeoff").apply(lambda df: (df.growth_rate > 1e-6).sum()).reset_index()

