#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 14:39:00 2022

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
    df=pd.read_csv(dec+'/carbon_sources.csv')
    df.set_index('carbon_sources',inplace=True)
    
    big=list(range(0,int(round(math.sqrt(len(df.index))))))
    lst = [list(z) for z in itertools.product(big, repeat=2)]
    sns.set(font_scale=2.2)
    fig, axes = plt.subplots(round(math.sqrt(len(df.index))), round(math.sqrt(len(df.index))), figsize=(23.5, 30))
    sns.set(rc={'axes.facecolor':'#f3f3f3', 'figure.facecolor':'white'})
    for i in range(len(df.index)):
        rows = df.loc[df.index[i]].T
        data=pd.DataFrame()
        data[rows.name]=[v for v in rows]
        sns.swarmplot(ax=axes[lst[i][0],lst[i][1]],data=data,palette = ['#e50808'],size=4)
        sns.boxplot(ax=axes[lst[i][0],lst[i][1]], data=data, palette=["#86c76a"],fliersize=2.5,showmeans=True,width=0.38)
        fig.savefig(dec+'/carbon_sources.png',bbox_inches='tight')
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    