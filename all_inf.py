#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 17:54:26 2022

@author: omidard
"""

from dgap import m9,exchange_reactions,essential_substrates, total_dir
import cobra
import pandas as pd
from cobra.io import load_json_model
from glob import glob
from cobra.manipulation.modify import rename_genes
from cobra import Model, Reaction, Metabolite
import os
import multiprocessing as mp

aux=[]
carbs=[]
exchanges = []
all_inf =pd.read_csv('/home/omidard/all_inf.csv')
for i in all_inf.columns[25:88]:
    if 'Unnamed' in i:
        all_inf = all_inf.drop(i,axis=1)
    if 'EX_' in i and 'eliminated_' not in i and 'growth_on_' not in i:
        exchanges.append(i)
    if 'growth_on_' in i:
        car = i.replace('growth_on_','')
        carbs.append(car)
    if 'eliminated_' in i:
        ax = i.replace('eliminated_','')
        aux.append(ax)



for i in range(len(all_inf)):
    model = load_json_model(all_inf['file'][i]+'.json')
    m9(model)
    fx = model.optimize().fluxes
    all_inf['wild_type_growth_rate'][i] = fx.BIOMASS2
    for e in exchanges:
        all_inf[e][i] = fx[e]
    for c in carbs:
        model = load_json_model(all_inf['file'][i]+'.json')
        m9(model)
        if 'EX_glc_D_e' not in c:
            model.reactions.EX_glc_D_e.lower_bound = 0
            model.reactions.get_by_id(c).lower_bound = 15
        model.reactions.get_by_id(c).lower_bound = 15
        fxc = model.optimize().fluxes
        all_inf['growth_on_'+c][i] = fxc.BIOMASS2
    for a in aux:
        model = load_json_model(all_inf['file'][i]+'.json')
        m9(model)
        model.reactions.get_by_id(e).lower_bound = 0
        fxa = model.optimize().fluxes
        all_inf['eliminated_'+a][i] = fxa.BIOMASS2
all_inf.to_csv('/home/omidard/all_informations.csv')
        
        
        
        
        
        
        
        
        
    
    
    
    