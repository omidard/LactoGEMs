#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  9 09:28:28 2022

@author: omidard
"""

from dgap import m9,gap_extract,none_ess_gap
import cobra
import pandas as pd
from cobra.io import load_json_model
from glob import glob
from cobra.manipulation.modify import rename_genes
from cobra import Model, Reaction, Metabolite
import os



names = []
rootdir = '/home/omidard/allgems'
for file in os.listdir(rootdir):
    local = os.path.join(rootdir, file)
    if os.path.isdir(local):
        names.append(local)

for name in names:
    os.makedirs(name+'/corr',exist_ok=True)
    directory1 = name+'/feasible/'
    directory2 = name+'/corr/'
    models =glob('%s/*.json'%directory1)
    for mod in models:
        model = load_json_model(mod)
        total_gaps = gap_extract(model)
        none_ess_gap(total_gaps,directory1,directory2)
