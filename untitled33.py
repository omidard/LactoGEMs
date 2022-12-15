#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 17:46:13 2022
test pc - plot is working
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

models =glob('%s/*.json'%'/Users/omidard/Desktop/lacba/allgems/allgems/allgems/inhouse/corr/')
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
        species_exchange_fluxes.to_csv('/Users/omidard/Desktop/lacba/allgems/allgems/allgems/inhouse/species_exchange_fluxes.csv')


#%%plot
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
        #axes[lst[i][0],lst[i][1]].set_facecolor('#c7c7c7')
        palette = ["#cc0000"]
    if len(con1)>0 and len(con2)==0:
        palette = ["#427fff"]
    if len(con1)==0 and len(con2)>0:
        palette = ["#4ea775"]
    sns.boxplot(ax=axes[lst[i][0],lst[i][1]], data=data, palette=palette,fliersize=2.5,showmeans=True,width=0.38)
    axes[lst[i][0],lst[i][1]].set_title(str(len(rows)-len(con3))+' strains', loc='left', y=0.33, x=0.15, fontsize=13,rotation=90)
    fig.savefig('/Users/omidard/Desktop/lacba/allgems/allgems/allgems/inhouse/pc2.png',bbox_inches='tight')  
    
#%%%

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

                         
#%%
dec = '/Users/omidard/Desktop/lacba/allgems/allgems/allgems/inhouse'
genera = dec.replace('/Users/omidard/Desktop/lacba/allgems/allgems/allgems/','')
df=pd.read_csv(dec+'/carbon_sources.csv')
df.set_index('carbon_sources',inplace=True)


#%%%
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
    y=0.15+i/15
    xy=(x,y)
    ax.annotate(species_wide_auxotrophies[i],xy=xy,xycoords='axes fraction')
palette=sns.color_palette("Set2")
sns.heatmap(df, annot=False, linewidths=.5, ax=ax,cmap="gist_gray_r")
fig.savefig(dec+'/'+genera+'auxotrophies',bbox_inches='tight')


#%% community test
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

data=pd.read_csv('/Users/omidard/Desktop/lacba/gems/all_inf2.csv')
model_db = '/Users/omidard/Desktop/lacba/gems/gems2'
manifest = build(data,out_folder="models",model_db=None,cutoff=0.0001,threads=60)
com = load_pickle("models/sample_1.pickle")
print(len(com.reactions))
medium= pd.read_csv('/Users/omidard/Desktop/lacba/gems/medium.csv')
res = grow(manifest, model_folder="models", medium=medium, tradeoff=0.5, threads=2)
tradeoff_rates = tradeoff(manifest, model_folder="models", medium=medium, threads=2)
tradeoff_rates.head()
tradeoff_rates.groupby("tradeoff").apply(lambda df: (df.growth_rate > 1e-6).sum()).reset_index()
























   
    
    
    
    
    