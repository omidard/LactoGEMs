#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 11:11:08 2022
find genes difference between best producers and non producers
@author: omidard
"""

from cobra.io import load_json_model,save_json_model
from dgap import m9
import cobra
import pandas as pd
from cobra.manipulation.modify import rename_genes
from cobra import Model, Reaction, Metabolite



exchanges = ['EX_ac_e','EX_actn_S_e','EX_akg_e',
'EX_for_e','EX_fum_e','EX_h2s_e','EX_mnl_e','EX_orn_e','EX_pyr_e','EX_succ_e',
'EX_lac_L_e','EX_13ppd_e','EX_4abut_e','EX_acald_e',
'EX_btd_RR_e','EX_etoh_e']
all_inf=pd.read_csv('/Users/omidard/Desktop/lacba/gems/all_informations.csv')

all_inf2=all_inf
bpc=[]
for i in exchanges:
    bp = all_inf2.sort_values(by=[i],ascending=False)
    df = pd.DataFrame()
    df['Species']=bp['Species'][:10]
    df['Strain']=bp['strain'][:10]
    df['id']=bp['id'][:10]
    df[i]=bp[i][:10]
    name = [' '.join(i) for i in zip(bp["Species"][:10].map(str),bp["strain"][:10])]
    df['name'] = name
    bpc.append(df)


all_inf2=all_inf
wpc=[]
for i in exchanges:
    wp = all_inf2.sort_values(by=[i],ascending=True)
    df = pd.DataFrame()
    df['Species']=wp['Species'][:10]
    df['Strain']=wp['strain'][:10]
    df['id']=wp['id'][:10]
    df[i]=wp[i][:10]
    name = [' '.join(i) for i in zip(wp["Species"][:10].map(str),wp["strain"][:10])]
    df['name'] = name
    wpc.append(df)
    
    
tgc=[]  #list of active genes in producers but not in non-producers  
for i in range(len(bpc)):
    nzfx = []
    for mo in bpc[i]['id']:
        mod = mo+'.json'
        model = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor/'+mod)
        m9(model)
        fx = model.optimize().fluxes
        for j in range(len(fx)):
            if fx[j] != 0:
                if fx.index[j] not in nzfx:
                    nzfx.append(fx.index[j])
    nzfx2=[]                
    for mo in wpc[i]['id']:
        mod = mo+'.json'
        model = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor/'+mod)
        m9(model)
        fx = model.optimize().fluxes
        for v in range(len(fx)):
            if fx[v] != 0:
                nzfx2.append(fx.index[v])
    
    targetgenes=[]
    for b in nzfx:
        if b not in nzfx2:
            targetgenes.append(b)
    df=pd.DataFrame()
    df[bpc[i].columns[3]] = targetgenes
    tgc.append(targetgenes)
    
#%% find effective genes on products out of tgc
for i in range(len(bpc)):
    for mo in bpc[i]['id']:
        mod = mo+'.json'
        model = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor/'+mod)
        m9(model)
        fx = model.optimize().fluxes
        modre=[]
        for re in model.reactions:
            modre.append(re.id)
        for j in tgc[i]:
            if j in modre:
                model.reactions.get_by_id(j).lower_bound = 0
                model.reactions.get_by_id(j).upper_bound = 0
                
            
        
    
    
    
    
    
    
    
    
    
    
    
    