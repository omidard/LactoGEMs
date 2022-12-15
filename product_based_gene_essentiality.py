#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 24 10:34:58 2022
inv gene analysis
@author: omidard
"""

from glob import glob
from cobra.io import load_json_model,save_json_model
from dgap import m9
import cobra
import pandas as pd
from cobra.manipulation.modify import rename_genes
from cobra import Model, Reaction, Metabolite

#%%
dfc = []
files =glob('%s/*.csv'%'/Users/omidard/Desktop/lacba/allinf/invsuc')
for i in files:
    df =pd.read_csv(i)
    dfc.append(df)
for i in dfc:
    i.set_index('gene',inplace=True)

inv_gen = pd.concat(dfc,axis=1)
inv_gen.fillna(0, inplace=True)
#inv_gen.to_csv('/Users/omidard/Desktop/lacba/allinf/inv/inv_gen.csv')
#%% strains
strains=[]
for i in inv_gen.columns:
    if 'wildtype' in i:
        st=i.replace('wildtype','')
        strains.append(st)
for st in strains:        
    for i in inv_gen.index:
        if inv_gen['knocked_out'+st][i] == '0.0':
            inv_gen['knocked_out'+st][i] = 0
            
for st in strains:        
    for i in inv_gen.index:
        if inv_gen['knocked_out'+st][i] != 'essential':
            inv_gen['knocked_out'+st][i] = int(float(inv_gen['knocked_out'+st][i]))
            
#%%
collect = []
for st in strains:
    pos = []
    posval = []
    neg = []
    negval=[]
    non = []
    nonval=[]
    for i in inv_gen.index:
        if inv_gen['wildtype'+st][i] != 0:
            if inv_gen['knocked_out'+st][i] != 'essential':
                if inv_gen['wildtype'+st][i] == inv_gen['knocked_out'+st][i]:
                    non.append(i)
                    nonval.append(0)
                if inv_gen['wildtype'+st][i] > inv_gen['knocked_out'+st][i]:
                    neg.append(i)
                    negval.append((inv_gen['knocked_out'+st][i] - inv_gen['wildtype'+st][i])/inv_gen['wildtype'+st][i]*100)
                if inv_gen['wildtype'+st][i] < inv_gen['knocked_out'+st][i]:
                    pos.append(i)
                    posval.append((inv_gen['knocked_out'+st][i] - inv_gen['wildtype'+st][i])/inv_gen['wildtype'+st][i]*100)
    dfpo = pd.DataFrame()
    dfpo['gene'] = pos
    dfpo['effect'] = posval
    dfpo['model'] = [st for i in range(len(dfpo))]
    dfpo.set_index('model',inplace=True)
    
    dfne = pd.DataFrame()
    dfne['gene']=neg
    dfne['effect']=negval
    dfne['model'] = [st for i in range(len(dfne))]
    dfne.set_index('model',inplace=True)
    
    dfno = pd.DataFrame()
    dfno['gene']=non
    dfno['effect']= nonval
    dfno['model'] = [st for i in range(len(dfno))]
    dfno.set_index('model',inplace=True)
    dft = pd.concat([dfpo,dfno,dfne],axis=0)
    collect.append(dft)
    
#%%
gg = pd.concat(collect,axis=0)
gg.sort_values(by='effect', ascending=True,inplace=True)
gg.to_csv('/Users/omidard/Desktop/lacba/allinf/invsuc/suc_gene_effects.csv')  
#%%% most repeated genes
gg['index'] = [v for v in range(len(gg))]
gg['model'] = [v for v in gg.index]
gg.set_index('index',inplace=True)
#%%
for i in gg.index:
    if gg['effect'][i] != -100:
        gg.drop(i, axis=0,inplace=True)
#%%
zz=['EX_glc_D_e','EX_h_e','GLCpts','EX_h2o_e']
tt=[]
ff = gg['gene'].value_counts()[:100]
#ff = gg['gene'].value_counts()
for i in ff.index:
    tt.append(i)
top10 = gg

for i in top10.index:
    if top10['gene'][i] not in tt:
        top10.drop(i, axis=0,inplace=True)
for i in top10.index:
    if top10['gene'][i] in zz:
        top10.drop(i, axis=0,inplace=True)
#%%

perc = pd.DataFrame()
perc['gene']=[i for i in ff.index]
perc['presence'] = [i/len(strains)*100 for i in ff]
#%%
top10['frequency'] = [0 for i in range(len(top10))]

for i in top10.index:
    for j in perc.index:
        if top10['gene'][i] == perc['gene'][j]:
            top10['frequency'][i] = perc['presence'][j]
#%%
all_inf=pd.read_csv('/Users/omidard/Desktop/lacba/gems/all_informations.csv')
top10['genus'] = [0 for i in range(len(top10))]
top10['species'] = [0 for i in range(len(top10))]
for i in top10.index:
    for j in all_inf.index:
        if top10['model'][i] == all_inf['id'][j]+'.json':
            top10['genus'][i] = all_inf['genus'][j]
            top10['species'][i] = all_inf['Species'][j]
#%%
top10.to_csv('/Users/omidard/Desktop/lacba/allinf/invsuc/suc_essential_genes.csv')

#%%
from dgap import m9
models =glob('%s/*.json'%'/Users/omidard/Desktop/biomassed/')
for mod in models:
    model = load_json_model(mod)
    m9(model)
    print(len(model.reactions),len(model.genes),model.optimize().fluxes.BIOMASS2)
    















    
        

