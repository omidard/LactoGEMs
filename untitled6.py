#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 19:26:32 2022

@author: omidard
"""



import cobra
import pandas as pd
from cobra.io import load_json_model
#%%
all_inf=pd.read_csv('/Users/omidard/Desktop/lacba/inf.csv')
models=[]
for i in all_inf['id'][:4]:
    mod=i+'.json'
    model=load_json_model('/Users/omidard/Desktop/lacba/gems/allgems/'+mod)
    models.append(model)
#%%
standardized_models =[]
for model in models:
    standardized_models.append(set_default_configs_and_snm3_medium(model))
#%%
for model in standardized_models:
    print(f"Growth {model.id}: {model.slim_optimize()}")
#%%
for model in standardized_models:
    old_medium = model.medium
    m2, extension = gapfill_medium(model)
    print(f"We need to add following metabolites to the medium: {extension}")
    print(f"Growth {m2.id}: {m2.slim_optimize()}")
    m2.medium = old_medium # Rest for next step