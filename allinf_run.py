#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 12 11:12:08 2022

@author: omidard
"""

from dgap import m9,allinf
import cobra
import pandas as pd
from cobra.io import load_json_model
from glob import glob
from cobra.manipulation.modify import rename_genes
from cobra import Model, Reaction, Metabolite
import os
import multiprocessing as mp
import pickle


models =glob('%s/*.json'%'/home/omidard/gems/allgems/')
pool = mp.Pool(mp.cpu_count())
dfc = pool.map(allinf, [modelid for modelid in models[0:240]])
all_informations = pd.concat(dfc,axis=0)
all_informations.to_csv('/home/omidard/all_informations.csv')
