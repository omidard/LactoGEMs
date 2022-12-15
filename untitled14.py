#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 20:05:44 2022

@author: omidard
"""
import cobra
from dgap import m9
from cobra.io import load_json_model
from glob import glob



directory = '/Users/omidard/Desktop/lacba/allgems/allgems/Weissella_confusa/biomassed'
models =glob('%s/*.json'%directory)

for mod in models:
    model=load_json_model(mod)
    m9(model)
    if model.optimize().fluxes.BIOMASS2 < 0.1:
        print(model.id,model.optimize().fluxes.BIOMASS2)
        cobra.io.json.save_json_model(model,'/Users/omidard/Desktop/lacba/allgems/allgems/Weissella_confusa/failed2/'+model.id)
    
    
    
#%%
import cobra
from dgap import m9
from cobra.io import load_json_model
from glob import glob

directory = '/Users/omidard/Desktop/lacba/allgems/allgems/Pediococcus_acidilactici/feasible'
models =glob('%s/*.json'%directory)

for mod in models:
    model=load_json_model(mod)
    m9(model)
    print(model.id, model.optimize().fluxes.BIOMASS2)
    
    
#%%


directory = '/Users/omidard/Desktop/lacba/allgems/allgems/Ligilactobacillus_salivarius/failed2'
directory2 = '/Users/omidard/Desktop/lacba/allgems/allgems/Ligilactobacillus_salivarius/feasible2'
models =glob('%s/*.json'%directory)
models2 =glob('%s/*.json'%directory2)
fis = []
for mod in models:
    model=load_json_model(mod)
    fis.append(model.id)
    
ttl = []
for mod in models2:
    model=load_json_model(mod)
    ttl.append(model.id)


for i in ttl:
    if i not in fis:
        print(i)

#%%

import cobra
import pandas as pd
import seaborn as sns
from cobra.io import load_json_model
from glob import glob
from cobra.manipulation.modify import rename_genes
from cobra import Model, Reaction, Metabolite
#model = load_json_model('/Users/omidard/Desktop/lacba/allgems/allgems/inhouse/gapfilled/PAL008.json.json')




reaction = Reaction('BIOMASS2')
reaction.name = 'Biomass production'
reaction.subsystem = 'Biomass'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default

#metabolites
thmpp_c = model.metabolites.get_by_id('thmpp_c')
CPS_LBR_c = model.metabolites.get_by_id('CPS_LBR_c')
DNA_LBR_c = model.metabolites.get_by_id('DNA_LBR_c')
LIP_LBR_c = model.metabolites.get_by_id('LIP_LBR_c')
LTAtotal_LBR_c = model.metabolites.get_by_id('LTAtotal_LBR_c')
PGlac2_c = model.metabolites.get_by_id('PGlac2_c')
PROT_LBR_c = model.metabolites.get_by_id('PROT_LBR_c')
RNA_LBR_c = model.metabolites.get_by_id('RNA_LBR_c')
atp_c = model.metabolites.get_by_id('atp_c')
btn_c = model.metabolites.get_by_id('btn_c')
coa_c = model.metabolites.get_by_id('coa_c')
h2o_c = model.metabolites.get_by_id('h2o_c')
nad_c = model.metabolites.get_by_id('nad_c')
pydx5p_c = model.metabolites.get_by_id('pydx5p_c')
thf_c = model.metabolites.get_by_id('thf_c')
udcpdp_c = model.metabolites.get_by_id('udcpdp_c')
adp_c = model.metabolites.get_by_id('adp_c')
h_c = model.metabolites.get_by_id('h_c')
pi_c = model.metabolites.get_by_id('pi_c')

reaction.add_metabolites({
    thmpp_c: -0.0001,
    CPS_LBR_c: -0.078,
    DNA_LBR_c: -0.205,
    atp_c: -41.2,
    LIP_LBR_c: -0.106,
    LTAtotal_LBR_c: -0.006,
    PGlac2_c: -0.009,
    PROT_LBR_c: -3.311,
    RNA_LBR_c: -0.926,
    btn_c: -0.00001,
    coa_c: -0.002,
    h2o_c: -41.2,
    nad_c: -0.002,
    pydx5p_c: -0.000001,
    thf_c: -0.00001,
    udcpdp_c: -0.002,
    adp_c: 41.2,
    h_c: 41.2,
    pi_c: 41.2
})
reaction.gene_reaction_rule = '(BIOMASS)'
reaction.reaction
model.add_reactions([reaction])
model.repair()

cobra.io.json.save_json_model(model,'/Users/omidard/Desktop/lacba/allgems/allgems/inhouse/failed2'+model.id)

#%%%

model = load_json_model('/Users/omidard/Desktop/lacba/allgems/allgems/inhouse/failed2/PAL008.json')
m9(model)
print(model.optimize().fluxes.BIOMASS2)


#%%

directory = '/Users/omidard/Desktop/lacba/allgems/allgems/Lacticaseibacillus_rhamnosus/biomassed'
models =glob('%s/*.json'%directory)

for mod in models:
    model=load_json_model(mod)
    cobra.io.json.save_json_model(model,'/Users/omidard/Desktop/lacba/allgems/allgems/Lacticaseibacillus_rhamnosus/feasible/'+model.id)
    




