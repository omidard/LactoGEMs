#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 18:16:29 2022
models debugging
@author: omidard
"""
# lactate transporter should be added ----done
# GLCpts and GLCabc from ref should be added to all models
# model.reactions.get_by_id('L-LACD').id='L_LACD' ---- done
# reaction.remove_from_model() EX_ptrc_e add EX_ptrc_e ---done
#%%
from cobra.io import load_json_model,save_json_model
from dgap import m9
import cobra
import pandas as pd
from cobra.manipulation.modify import rename_genes
from cobra import Model, Reaction, Metabolite
all_inf=pd.read_csv('/Users/omidard/Desktop/lacba/inf.csv')

#%%
target=[]
targetid=[]
for mo in all_inf['id']:
    mod = mo+'.json'
    model = load_json_model('/Users/omidard/Desktop/lacba/gems/allgems/'+mod)
    for re in model.reactions:
        if 'LACLt' in re.id or 'LACDt_c' in re.id:
            print(re.id)
            target.append(re)
            targetid.append(re.id)
            
            
#%% general reactions
reactions=[]
ref = load_json_model('/Users/omidard/Desktop/lacba/ref_model/COBRAModel.json')
miss = ['GAP','Diffusion','SINK','DEMAND','EXCHANGE','ORPHAN','spontaneous']
for i in miss:
    re = ref.genes.get_by_id(i).reactions
    for v in re:
        reactions.append(v)
        
for mo in all_inf['id']:
    mod = mo+'.json'
    model = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor/'+mod)
    model.add_reactions(reactions)
    save_json_model(model,'/Users/omidard/Desktop/lacba/gems/allcor/'+model.id)
#%% GLCpts and GLCabc from ref should be added to all models
mreactions=[]
ref = load_json_model('/Users/omidard/Desktop/lacba/ref_model/COBRAModel.json')
miss = ['GLCpts','GLCabc']
for i in miss:
    reaction = ref.reactions.get_by_id(i)
    reaction2 = Reaction(reaction.id)
    reaction2.name = reaction.name
    reaction2.subsystem = reaction.subsystem
    reaction2.lower_bound = reaction.lower_bound
    reaction2.upper_bound = reaction.upper_bound
    reaction2.add_metabolites(reaction.metabolites)
    reaction2.gene_reaction_rule = '(GAP)'
    mreactions.append(reaction2)

for mo in all_inf['id']:
    mod = mo+'.json'
    model = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor/'+mod)
    model.add_reactions(mreactions)
    save_json_model(model,'/Users/omidard/Desktop/lacba/gems/allcor/'+model.id)


#%% correction of EX_ptrc_e directionality
for mo in all_inf['id']:
    mod = mo+'.json'
    model = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor/'+mod)
    model.reactions.EX_ptrc_e.remove_from_model()
    model.add_boundary(model.metabolites.get_by_id("ptrc_e"), type="exchange")
    print('added to ',' ',model.id)
    save_json_model(model,'/Users/omidard/Desktop/lacba/gems/allcor/'+model.id)















