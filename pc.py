#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 12 23:03:12 2022
production-consumption
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


#'/home/omidard/allgems'
dirs = total_dir('/Users/omidard/Desktop/lacba/allgems/allgems/allgems')
for dec in dirs:
    directory = dec+'/corr'
    print(directory)
    models =glob('%s/*.json'%directory)
    flux_collector=[]
    for mod in models:
        model = load_json_model(mod)
        try:
            flux_frame = exchanges_fluxes(model)
            flux_collector.append(flux_frame)
        except:
            pass
    species_exchange_fluxes = pd.concat(flux_collector, axis=1)
    for i in species_exchange_fluxes.index:
        rows = species_exchange_fluxes.loc[i]
        if all([ v == 0 for v in rows ]):
            species_exchange_fluxes = species_exchange_fluxes.drop(index=(i))
            species_exchange_fluxes.to_csv(dec+'/species_exchange_fluxes.csv')
    
    
    #plot
    big=list(range(0,int(round(math.sqrt(len(species_exchange_fluxes.index))))+1))
    lst = [list(z) for z in itertools.product(big, repeat=2)]
    fig, axes = plt.subplots(round(math.sqrt(len(species_exchange_fluxes.index)))+1, round(math.sqrt(len(species_exchange_fluxes.index)))+1, figsize=(23.5, 30))
    sns.set(rc={'axes.facecolor':'#f3f3f3', 'figure.facecolor':'white'})
    for i in range(len(species_exchange_fluxes.index)):
        rows = species_exchange_fluxes.loc[species_exchange_fluxes.index[i]].T
        data=pd.DataFrame()
        data[rows.name]=[v for v in rows]
        
        con1=[]
        con2=[]
        con3=[]
        for v in rows:
            if v>0:
                con1.append(1)
            if v<0:
                con2.append(1)
            if v ==0:
                con3.append(1)
                
        sns.swarmplot(ax=axes[lst[i][0],lst[i][1]],data=data,palette = ['#000000'],size=1.5)
        if len(con1)>0 and len(con2)>0:
            palette = ["#cc0000"]
        if len(con1)>0 and len(con2)==0:
            palette = ["#427fff"]
        if len(con1)==0 and len(con2)>0:
            palette = ["#4ea775"]
        sns.boxplot(ax=axes[lst[i][0],lst[i][1]], data=data, palette=palette,fliersize=2.5,showmeans=True,width=0.38)
        axes[lst[i][0],lst[i][1]].set_title(str(len(rows)-len(con3))+' strains', loc='left', y=0.33, x=0.15, fontsize=13,rotation=90)
        fig.savefig(dec+'/pc.png',bbox_inches='tight')








