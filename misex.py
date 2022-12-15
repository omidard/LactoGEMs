#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 13 12:27:09 2022
add missing exchanges
@author: omidard
"""

from dgap import m9,exchanges_fluxes,total_dir
import cobra
import pandas as pd
from cobra.io import load_json_model
from glob import glob
from cobra.manipulation.modify import rename_genes
from cobra import Model, Reaction, Metabolite
import os
import matplotlib.pyplot as plt
import seaborn as sns
import itertools
import math

ref_dir = '/home/omidard/ref_model_dir/COBRAModel.json'
ref = load_json_model(ref_dir)
ref_ex=[]
for re in ref.reactions:
    if 'EX_' in re.id:
        ref_ex.append(re.id)
        
dirs = total_dir('/home/omidard/allgems')
for dec in dirs:
    directory = dec+'/corr'
    models =glob('%s/*.json'%directory)
    mod_ex=[]
    for mod in models:
        model = load_json_model(mod)
        for re in model.reactions:
            if 'EX_' in re.id:
                mod_ex.append(re.id)
        missing=[]                
        for r in ref_ex:
            if r not in mod_ex:
                missing.append(r)
        reactions=[]
        for re in missing:
            reactions.append(ref.reactions.get_by_id(re))
        model.add_reactions(reactions)
        model.repair()
        cobra.io.json.save_json_model(model,directory+'/'+model.id)
        
        
                