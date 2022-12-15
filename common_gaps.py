#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 19 12:28:55 2022

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





all_gaps_count=pd.read_csv('/Users/omidard/Desktop/lacba/gems/all_gaps_count.csv')
all_gaps_count.set_index('gaps',inplace=True)
#%% turn it into percentage table
for col in all_gaps_count.columns:
    column = all_gaps_count[col]
    max_value = column.max()
    for i in range(len(column)):
        all_gaps_count[col][i]=all_gaps_count[col][i]/max_value*100
#%%% fill nan with zero
import numpy as np
df = all_gaps_count.fillna(0)

#%%%
df.to_csv('/Users/omidard/Desktop/lacba/gems/all_gaps_count2.csv')
#%%
df2=df.sort_values(by = [c for c in df.columns],ascending = False)
#%%


import matplotlib.pyplot as plt
import seaborn as sns
sns.set_theme()

# Load the example flights dataset and convert to long-form
#flights_long = sns.load_dataset("flights")
#flights = flights_long.pivot("month", "year", "passengers")

# Draw a heatmap with the numeric values in each cell
plt.xlabel('gaps', fontsize=19)
plt.xlabel('gaps', fontsize=16)
f, ax = plt.subplots(figsize=(30, 15))
plt.xlabel('gaps', fontsize=16)
sns.heatmap(df2[:50].T, linewidths=.5, ax=ax,robust=True)
f.savefig('/Users/omidard/Desktop/lacba/gems/allgaps2.png',bbox_inches='tight')
#%%

import pandas as pd
import matplotlib.pyplot as plt
from catscatter import catscatter
dfc=[]
for i in range(40):
    rows = df2.loc[df.index[i]]
    df3=pd.DataFrame()
    df3['gap_count']=rows
    df3['gaps']=[rows.name for i in range(len(rows))]
    dfc.append(df3)
df4 = pd.concat(dfc, axis=0)
#%%


df4['genera']=[i for i in df4.index]






kcolors=['#F73972','#F2B3C6','#144962']
# create the plot
plt.figure(figsize=(80,40))
catscatter(df4,'gaps','genera','gap_count',color=kcolors, ratio=20)
plt.xticks(fontsize=35)
plt.yticks(fontsize=35)
plt.show()

































