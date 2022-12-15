#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 08:34:31 2022
auxotrophies plot
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



dirs = total_dir('/home/omidard/allgems')
for dec in dirs:
    genera = dec.replace('/home/omidard/allgems/','')
    df=pd.read_csv(dec+'/minimal_media.csv')
    df.set_index('metabolites',inplace=True)
    #added
    mostc=[]
    for i in range(len(df.index)):
        rows = df.loc[df.index[i]]
        most=[]
        for v in rows:
            if v==0:
                most.append(1)
        mostc.append(most)
    species_wide_auxotrophies=[]
    for i in range(len(mostc)):
        if len(mostc[i]) == len(df.columns):
            species_wide_auxotrophies.append(df.index[i])

    
    sns.set(font_scale=3)
    fig, ax = plt.subplots(figsize=(80,35))
    #added
    ax.annotate('species wide auxotrophies',xy=(1.15,0.85),xycoords='axes fraction')
    for i in range(len(species_wide_auxotrophies)):
        x=1.15
        y=0.15+i/17
        xy=(x,y)
        ax.annotate(species_wide_auxotrophies[i],xy=xy,xycoords='axes fraction')
    #end
    sns.heatmap(df, annot=False, linewidths=.5, ax=ax)
    fig.savefig(dec+'/'+genera+'auxotrophies',bbox_inches='tight')







