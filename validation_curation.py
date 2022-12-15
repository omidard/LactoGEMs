#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 22 04:11:06 2022

@author: omidard
"""

from cobra.io import load_json_model,save_json_model
from dgap import m9
import cobra
import pandas as pd
from cobra.manipulation.modify import rename_genes
from cobra import Model, Reaction, Metabolite



sod=pd.read_excel('/Users/omidard/Desktop/lacba/gems/single_omission.xlsx')
all_inf=pd.read_csv('/Users/omidard/Desktop/lacba/gems/all_informations.csv')
#%% plantarum WCFS1
for i in range(len(all_inf)):
    if 'LC2W' in all_inf['strain'][i]:
        mod = all_inf['id'][i]+'.json'
        model = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor/'+mod)
        m9(model)
        orgin = model.optimize().fluxes.BIOMASS2
        orgi = orgin/10
        org = orgi * 4
        df = pd.DataFrame()
        for j in sod['component'][0:50]:
            model = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor/'+mod)
            m9(model)
            model.reactions.get_by_id(j).lower_bound = 0
            bm = model.optimize().fluxes.BIOMASS2
            if bm < org:
                gr = 0
            if bm >= org:
                gr = 1
            df[j]=[gr]

df.index = ['Prediction']
#%%

df2 = pd.DataFrame()
for i in range(len(sod['casei_LC2W'][0:50])):
    if sod['casei_LC2W'][i] == 'N':
        gr = 1
    if sod['casei_LC2W'][i] == 'E':
        gr = 0
    if sod['casei_LC2W'][i] == 'ND':
        gr = 1
    df2[sod['component'][i]] = [gr]
df2.index = ['Experimental']
#%%
df3 = pd.concat([df,df2],axis=0)

#%% models curation based on results

for i in range(len(all_inf)):
    if 'plantarum' in all_inf['Species'][i]:
        mod = all_inf['id'][i]+'.json'
        model = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor/'+mod)
        
        biomass = model.reactions.get_by_id("BIOMASS2")
        
        
        #remove metabolites
        biomass.subtract_metabolites({model.metabolites.get_by_id("pydx5p_c"): -0.000001})
        biomass.subtract_metabolites({model.metabolites.get_by_id("btn_c"): -0.00001})
        biomass.subtract_metabolites({model.metabolites.get_by_id("thmpp_c"): -0.0001})
        #add metabolites
        biomass.add_metabolites({model.metabolites.get_by_id("ascb6p_c"): -0.00001})
        biomass.add_metabolites({model.metabolites.get_by_id("ribflv_c"): -0.00001})
        model.repair()
        
        cobra.io.json.save_json_model(model,'/Users/omidard/Desktop/lacba/gems/allcor2/'+model.id)
        
#%%

for i in range(len(all_inf)):
    if 'mesenteroides' in all_inf['Species'][i]:
        mod = all_inf['id'][i]+'.json'
        model = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor/'+mod)
        biomass = model.reactions.get_by_id("BIOMASS2")
        #remove metabolites
        #biomass.subtract_metabolites({model.metabolites.get_by_id("coa_c"): -0.002})
        biomass.subtract_metabolites({model.metabolites.get_by_id("btn_c"): -0.00001})
        biomass.subtract_metabolites({model.metabolites.get_by_id("thmpp_c"): -0.0001})
        #add metabolites
        #biomass.add_metabolites({model.metabolites.get_by_id("ascb6p_c"): -0.00001})
        #biomass.add_metabolites({model.metabolites.get_by_id("ribflv_c"): -0.00001})
        model.repair()
        cobra.io.json.save_json_model(model,'/Users/omidard/Desktop/lacba/gems/allcor2/'+model.id)         
#%%
for i in range(len(all_inf)):
    if 'delbrueckii' in all_inf['Species'][i]:
        mod = all_inf['id'][i]+'.json'
        model = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor/'+mod)
        biomass = model.reactions.get_by_id("BIOMASS2")
        #remove metabolites
        biomass.subtract_metabolites({model.metabolites.get_by_id("pydx5p_c"): -0.000001})
        biomass.subtract_metabolites({model.metabolites.get_by_id("btn_c"): -0.00001})
        biomass.subtract_metabolites({model.metabolites.get_by_id("thmpp_c"): -0.0001})
        #add metabolites
        biomass.add_metabolites({model.metabolites.get_by_id("pydxn_c"): -0.00001})
        #biomass.add_metabolites({model.metabolites.get_by_id("ribflv_c"): -0.00001})
        model.repair()
        cobra.io.json.save_json_model(model,'/Users/omidard/Desktop/lacba/gems/allcor2/'+model.id)         
#%%
for i in range(len(all_inf)):
    if 'salivarius' in all_inf['Species'][i]:
        mod = all_inf['id'][i]+'.json'
        model = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor/'+mod)
        biomass = model.reactions.get_by_id("BIOMASS2")
        #remove metabolites
        biomass.subtract_metabolites({model.metabolites.get_by_id("thf_c"): -0.00001})
        #biomass.subtract_metabolites({model.metabolites.get_by_id("pydx5p_c"): -0.000001})
        biomass.subtract_metabolites({model.metabolites.get_by_id("btn_c"): -0.00001})
        biomass.subtract_metabolites({model.metabolites.get_by_id("thmpp_c"): -0.0001})
        #add metabolites
        #biomass.add_metabolites({model.metabolites.get_by_id("pydxn_c"): -0.00001})
        biomass.add_metabolites({model.metabolites.get_by_id("ribflv_c"): -0.00001})
        model.repair()
        cobra.io.json.save_json_model(model,'/Users/omidard/Desktop/lacba/gems/allcor2/'+model.id)
#%%
for i in range(len(all_inf)):
    if 'paracasei' in all_inf['Species'][i]:
        mod = all_inf['id'][i]+'.json'
        model = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor/'+mod)
        biomass = model.reactions.get_by_id("BIOMASS2")
        #remove metabolites
        biomass.subtract_metabolites({model.metabolites.get_by_id("thf_c"): -0.00001})
        #biomass.subtract_metabolites({model.metabolites.get_by_id("pydx5p_c"): -0.000001})
        biomass.subtract_metabolites({model.metabolites.get_by_id("btn_c"): -0.00001})
        biomass.subtract_metabolites({model.metabolites.get_by_id("thmpp_c"): -0.0001})
        #add metabolites
        #biomass.add_metabolites({model.metabolites.get_by_id("pydxn_c"): -0.00001})
        #biomass.add_metabolites({model.metabolites.get_by_id("ribflv_c"): -0.00001})
        model.repair()
        cobra.io.json.save_json_model(model,'/Users/omidard/Desktop/lacba/gems/allcor2/'+model.id)
#%%
for i in range(len(all_inf)):
    if 'rhamnosus' in all_inf['Species'][i]:
        mod = all_inf['id'][i]+'.json'
        model = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor/'+mod)
        biomass = model.reactions.get_by_id("BIOMASS2")
        #remove metabolites
        biomass.subtract_metabolites({model.metabolites.get_by_id("thf_c"): -0.00001})
        #biomass.subtract_metabolites({model.metabolites.get_by_id("pydx5p_c"): -0.000001})
        biomass.subtract_metabolites({model.metabolites.get_by_id("btn_c"): -0.00001})
        biomass.subtract_metabolites({model.metabolites.get_by_id("thmpp_c"): -0.0001})
        #add metabolites
        biomass.add_metabolites({model.metabolites.get_by_id("pydxn_c"): -0.00001})
        #biomass.add_metabolites({model.metabolites.get_by_id("ribflv_c"): -0.00001})
        model.repair()
        cobra.io.json.save_json_model(model,'/Users/omidard/Desktop/lacba/gems/allcor2/'+model.id)
#%%
for i in range(len(all_inf)):
    if 'acidophilus' in all_inf['Species'][i]:
        mod = all_inf['id'][i]+'.json'
        model = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor/'+mod)
        biomass = model.reactions.get_by_id("BIOMASS2")
        #remove metabolites
        biomass.subtract_metabolites({model.metabolites.get_by_id("thf_c"): -0.00001})
        #biomass.subtract_metabolites({model.metabolites.get_by_id("pydx5p_c"): -0.000001})
        biomass.subtract_metabolites({model.metabolites.get_by_id("btn_c"): -0.00001})
        biomass.subtract_metabolites({model.metabolites.get_by_id("thmpp_c"): -0.0001})
        #add metabolites
        #biomass.add_metabolites({model.metabolites.get_by_id("pydxn_c"): -0.00001})
        #biomass.add_metabolites({model.metabolites.get_by_id("ribflv_c"): -0.00001})
        model.repair()
        cobra.io.json.save_json_model(model,'/Users/omidard/Desktop/lacba/gems/allcor2/'+model.id)
#%%
for i in range(len(all_inf)):
    if 'reuteri' in all_inf['Species'][i]:
        mod = all_inf['id'][i]+'.json'
        model = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor/'+mod)
        biomass = model.reactions.get_by_id("BIOMASS2")
        #remove metabolites
        #biomass.subtract_metabolites({model.metabolites.get_by_id("thf_c"): -0.00001})
        #biomass.subtract_metabolites({model.metabolites.get_by_id("pydx5p_c"): -0.000001})
        biomass.subtract_metabolites({model.metabolites.get_by_id("btn_c"): -0.00001})
        biomass.subtract_metabolites({model.metabolites.get_by_id("thmpp_c"): -0.0001})
        #add metabolites
        #biomass.add_metabolites({model.metabolites.get_by_id("pydxn_c"): -0.00001})
        #biomass.add_metabolites({model.metabolites.get_by_id("ribflv_c"): -0.00001})
        model.repair()
        cobra.io.json.save_json_model(model,'/Users/omidard/Desktop/lacba/gems/allcor2/'+model.id)
#%%
for i in range(len(all_inf)):
    if 'oeni' in all_inf['Species'][i]:
        mod = all_inf['id'][i]+'.json'
        model = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor/'+mod)
        biomass = model.reactions.get_by_id("BIOMASS2")
        #remove metabolites
        #biomass.subtract_metabolites({model.metabolites.get_by_id("thf_c"): -0.00001})
        #biomass.subtract_metabolites({model.metabolites.get_by_id("pydx5p_c"): -0.000001})
        biomass.subtract_metabolites({model.metabolites.get_by_id("btn_c"): -0.00001})
        biomass.subtract_metabolites({model.metabolites.get_by_id("thmpp_c"): -0.0001})
        #add metabolites
        #biomass.add_metabolites({model.metabolites.get_by_id("pydxn_c"): -0.00001})
        #biomass.add_metabolites({model.metabolites.get_by_id("ribflv_c"): -0.00001})
        model.repair()
        cobra.io.json.save_json_model(model,'/Users/omidard/Desktop/lacba/gems/allcor2/'+model.id)


#%%test
reactions = []
for re in model.genes.GAP.reactions:
    reactions.append(re.id)
for re in model.genes.ORPHAN.reactions:
    reactions.append(re.id)
#%%
short_list=[]
model = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor2/GCF_000014445.1.json')
for re in reactions:
    model = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor2/GCF_000014445.1.json')
    m9(model)
    model.reactions.get_by_id(re).lower_bound = 0
    model.reactions.get_by_id(re).upper_bound = 0
    if model.optimize().fluxes.BIOMASS2 > 0.1:
        print(model.optimize().fluxes.BIOMASS2,re)
        short_list.append(re)
        #reaction = model.reactions.get_by_id(re)
        #reaction.remove_from_model()
#%%

model = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor2/GCF_000014445.1.json')
for re in short_list:
    if 'ADCS' not in re:
        m9(model)
        up = model.reactions.get_by_id(re).upper_bound
        low =  model.reactions.get_by_id(re).lower_bound
        model.reactions.get_by_id(re).lower_bound = 0
        model.reactions.get_by_id(re).upper_bound = 0
        if model.optimize().fluxes.BIOMASS2 > 0.1:
            print(model.optimize().fluxes.BIOMASS2,re)
            reaction = model.reactions.get_by_id(re)
            reaction.remove_from_model()
        if model.optimize().fluxes.BIOMASS2 < 0.1:
            model.reactions.get_by_id(re).lower_bound = low
            model.reactions.get_by_id(re).upper_bound = up

#%%% aminoacids auxotrophy correction automatically 1- detect auxotrophies with fn prediction
#2- find models with correct auxotrophies for those amino acids (tn)

fn=[]
for i in df3.columns:
    if df3[i]['Experimental'] == 1:
        if df3[i]['Prediction'] == 0:
            fn.append(i)
df4c=[]
for i in range(len(all_inf)):
    if 'plantarum' in all_inf['Species'][i]:
        mod = all_inf['id'][i]+'.json'
        model = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor/'+mod)
        m9(model)
        gr = model.optimize().fluxes.BIOMASS2
        grs = gr/10
        grs2 = grs*4
        df4 = pd.DataFrame()
        for j in fn:
            models = []
            model = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor/'+mod)
            m9(model)
            model.reactions.get_by_id(j).lower_bound=0
            if model.optimize().fluxes.BIOMASS2 >= grs2:
                models.append(mod)
            if model.optimize().fluxes.BIOMASS2 < grs2:
                models.append('not found')
            df4[j] = models
            df4c.append(df4)
#%%
for i in df4c:
    for j in i.columns:
        if 'not found' not in i[j][0]:
            print(i[j][0],' = ', j)
                
#%% map differences between target model and template model
target = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor/GCF_000203855.3.json')
template = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor/GCF_000175355.1.json')

#reactions of target model
targre = []
for re in target.reactions:
    targre.append(re.id)
#reactions of template model
tempre = []
for re in template.reactions:
    tempre.append(re.id)
#reactions which are not in template but are in target
exces=[] # exessive reaction
for i in targre:
    if i not in tempre:
        exces.append(i)
#orphan and gap reactions of target model
gaporf=[] #all reactions belong to gap and orphan (allowed to be removed)
for re in target.genes.GAP.reactions:
    gaporf.append(re.id)
for re in target.genes.ORPHAN.reactions:
    gaporf.append(re.id)
#reactions which are not in template but in target and are orphan or gap = final best target to remove
exces2 = [] #exessive reaction belongs to gap and orphan reactions (allowed to be removed)
for x in exces:
    if x in gaporf:
        exces2.append(x)
        

        
#%% test if removing selected reactions could meet the auxotrophy condition
memory = pd.DataFrame()
growth = []
for i in gaporf:
    target = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor/GCF_000026505.1.json')
    m9(target)
    target.reactions.get_by_id(i).lower_bound = 0
    target.reactions.get_by_id(i).upper_bound = 0
    gr = target.optimize().fluxes.BIOMASS2
    grs = gr/10
    grs2 = grs*4
    print(gr)
    growth.append(gr)
memory['growth']=growth
memory['reaction'] = gaporf
#%%
newtarg = []
for i in range(len(memory)):
    if memory['growth'][i] >= 0.3:
        newtarg.append(memory['reaction'][i])
#%%#amino acid
for i in newtarg:
    target = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor/GCF_000026505.1.json')
    m9(target)
    target.reactions.get_by_id(i).lower_bound = 0
    target.reactions.get_by_id(i).upper_bound = 0
    target.reactions.EX_thr_L_e.lower_bound = 0
    if target.optimize().fluxes.BIOMASS2 <=0.128:
        print(i,' ', target.optimize().fluxes.BIOMASS2)
#%%
        target.remove_reactions([target.reactions.get_by_id(i)])
        target.repair()
        cobra.io.json.save_json_model(target,'/Users/omidard/Desktop/lacba/gems/allcor2/'+target.id)
#%%
found = ['ACONT','ICDHy'] #glu fn solution
target = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor/GCF_000203855.3.json')
for i in found:
    target.remove_reactions([target.reactions.get_by_id(i)])
    target.repair()
    m9(target)
    print(target.optimize().fluxes.BIOMASS2)
    cobra.io.json.save_json_model(target,'/Users/omidard/Desktop/lacba/gems/allcor2/'+target.id)




#%% cystein fn curation template delbrueckii crl581
# 1- find a template
mods = []
grate = []
for i in range(len(all_inf)):
    mod = all_inf['id'][i]+'.json'
    model = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor/'+mod)
    m9(model)
    model.reactions.get_by_id('EX_cys_L_e').lower_bound = 0
    print(mod,' ', model.optimize().fluxes.BIOMASS2)
    if model.optimize().fluxes.BIOMASS2 > 0:
        mods.append(mod)
        grate.append(model.optimize().fluxes.BIOMASS2)
tempfind = pd.DataFrame()
tempfind['id'] = mods
tempfind['growth'] = grate


#%%
tempre=[]
for i in range(len(all_inf)):
    if 'CRL581' in all_inf['strain'][i]:
        mod1 = all_inf['id'][i]+'.json'
        template = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor/'+mod1)
        for re in template.reactions:
            tempre.append(re.id)
modre=[]
for i in range(len(all_inf)):
    if 'WCFS1' in all_inf['strain'][i]:
        mod2 = all_inf['id'][i]+'.json'
        model = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor/'+mod2)
        for re in model.reactions:
            modre.append(re.id)

misre=[]
for j in tempre:
    if j not in modre:
        misre.append(j)

misdf=pd.DataFrame()
grr=[]
for r in misre:
    template = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor/'+mod1)
    m9(template)
    template.reactions.get_by_id(r).lower_bound = 0
    template.reactions.get_by_id(r).upper_bound = 0
    print(r,' ', template.optimize().fluxes.BIOMASS2)
    grr.append(template.optimize().fluxes.BIOMASS2)
misdf['reaction']=misre
misdf['growth']=grr

targmis=[]
for m in range(len(misdf['growth'])):
    if misdf['growth'][m] > 0.3:
        targmis.append(misdf['reaction'][m])
#%%
for t in targmis:
    template = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor/'+mod1)
    m9(template)
    template.reactions.get_by_id(r).lower_bound = 0
    template.reactions.get_by_id(r).upper_bound = 0
    template.reactions.get_by_id('EX_cys_L_e').lower_bound = 0
    if template.optimize().fluxes.BIOMASS2 <= 0.144:
        print(t,' ',template.optimize().fluxes.BIOMASS2)
    
#%%    
template = load_json_model('/Users/omidard/Desktop/lacba/ref_model/COBRAModel.json')

rec = ['SADT','SADT2','AMPSDH','SO3R','SULR'] #cys fn solution
reaction = template.reactions.get_by_id('AMPSDH')
reaction2 = Reaction(reaction.id)
reaction2.name = reaction.name
reaction2.subsystem = reaction.subsystem
reaction2.lower_bound = reaction.lower_bound
reaction2.upper_bound = reaction.upper_bound
reaction2.add_metabolites(reaction.metabolites)
reaction2.gene_reaction_rule = '(GAP)'
model.add_reactions([reaction2])

#%%TEST FOR SULFUR ASSIMILATION REACTIONS WITHIN MODELS

for i in range(len(all_inf)):
    if 'LC2W' in all_inf['strain'][i]:
        mod = all_inf['id'][i]+'.json'
        model = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor/'+mod)


template = load_json_model('/Users/omidard/Desktop/lacba/ref_model/COBRAModel.json')
target_reactions = []
# reserved reactionsto be watched 'ADSK','TRPAS1'
#rec = ['SULabc','SADT','SADT2','AMPSDH','SO3R','SULR'] #cysteine
#rec = ['CS'] #citrate ACONTa is reserved

rec = ['HSERTA','METB1'] #methionine -> also  add reaction methionine synthase 
for i in rec:
    reaction = template.reactions.get_by_id(i)
    reaction2 = Reaction(reaction.id)
    reaction2.name = reaction.name
    reaction2.subsystem = reaction.subsystem
    reaction2.lower_bound = reaction.lower_bound
    reaction2.upper_bound = reaction.upper_bound
    reaction2.add_metabolites(reaction.metabolites)
    reaction2.gene_reaction_rule = '(GAP)'
    target_reactions.append(reaction2)
model.add_reactions(target_reactions)
model.repair()
#cobra.io.json.save_json_model(model,'/Users/omidard/Desktop/lacba/gems/allcor2/'+model.id)
#%%
m9(model)
print(model.optimize())
m9(model)
model.reactions.EX_cys_L_e.lower_bound = 0
print(model.optimize())
#%%
mm = []
for re in model.metabolites.met_L_c.reactions:
    print(re.id,' ',model.reactions.get_by_id(re.id).reaction,' ',model.reactions.get_by_id(re.id).lower_bound,' ',model.reactions.get_by_id(re.id).upper_bound)
    mm.append(re.id)
#%%
tm =[]
for re in template.metabolites.met_L_c.reactions:
    print(re.id,' ',template.reactions.get_by_id(re.id).reaction,' ',template.reactions.get_by_id(re.id).lower_bound,' ',template.reactions.get_by_id(re.id).upper_bound)
    tm.append(re.id)
#%%

for i in range(len(all_inf)):
    if 'CRL581' in all_inf['strain'][i]:
        mod = all_inf['id'][i]+'.json'
        model1 = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor/'+mod)

for i in range(len(all_inf)):
    if 'WCFS1' in all_inf['strain'][i]:
        mod = all_inf['id'][i]+'.json'
        model2 = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor/'+mod)

m2re=[]
for re in model2.reactions:
    m2re.append(re.id)
m1re=[]
for re in model1.reactions:
    m1re.append(re.id)
m3re=[]    
for i in m2re:
    if i not in m1re:
        m3re.append(i)

tre=[]
for i in m3re:
    model2 = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor/'+mod)
    m9(model2)
    model2.reactions.get_by_id(i).lower_bound = 0
    model2.reactions.get_by_id(i).upper_bound = 0
    if model2.optimize().fluxes.BIOMASS2 != 0:
        model2.reactions.get_by_id('EX_cys_L_e').lower_bound = 0
        if model2.optimize().fluxes.BIOMASS2 == 0:
            tre.append(i)
            
#%%
m9(model)
m9(model2)
for i in m2re:
    if i in m1re:
        if model2.reactions.get_by_id(i).lower_bound != model.reactions.get_by_id(i).lower_bound:
            print('low', i,'lcw= ',model2.reactions.get_by_id(i).lower_bound,'crl= ',model.reactions.get_by_id(i).lower_bound)
        if model2.reactions.get_by_id(i).upper_bound != model.reactions.get_by_id(i).upper_bound:
            print('up', i,'lcw= ',model2.reactions.get_by_id(i).upper_bound,'crl= ',model.reactions.get_by_id(i).upper_bound)
    
    
#%%% met auxotrophy solution

t = load_json_model('/Users/omidard/Desktop/lacba/ref_model/COBRAModel.json')
#%%

for re in t.metabolites.get_by_id('ala_B_c').reactions:
    print(re.id,' = ',t.reactions.get_by_id(re.id).reaction)
    
#%%
for re in model.metabolites.get_by_id('ala_B_c').reactions:
    print(re.id,' = ',model.reactions.get_by_id(re.id).reaction)
    
#%% add a new reaction <Methionine synthase II (cobalamin-independent)>
###  5mthf_c + hcys_L_c → h_c + met_L_c + thf_c  
### suchms_c  cyst_L_c  rhcys_c
from cobra import Model, Reaction, Metabolite   
    
reaction = Reaction('METS')
reaction.name = 'Methionine synthase II (cobalamin-independent)'
reaction.subsystem = 'Cysteine and methionine metabolism'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default



thf_c = model.metabolites.get_by_id('thf_c')


reaction.add_metabolites({
    model.metabolites.get_by_id('5mthf_c'): -1,
    model.metabolites.get_by_id('hcys_L_c'): -1,
    model.metabolites.get_by_id('h_c'): 1,
    model.metabolites.get_by_id('met_L_c'): 1,
    model.metabolites.get_by_id('thf_c'): 1,
})
reaction.gene_reaction_rule = '(GAP)'
reaction.reaction
model.add_reactions([reaction])
model.repair()



    
#%% pnto_R_c tracking pnto_R_c =====> coa_c
"""
PNTK  =  atp_c + pnto_R_c --> 4ppan_c + adp_c + h_c
PPNCL2  =  4ppan_c + ctp_c + cys_L_c --> 4ppcys_c + cmp_c + h_c + ppi_c
PPCDC  =  4ppcys_c + h_c --> co2_c + pan4p_c    
PTPATi  =  atp_c + h_c + pan4p_c <=> dpcoa_c + ppi_c    
DPCOAK  =  atp_c + dpcoa_c --> adp_c + coa_c + h_c    

#two strategies 1- pnto_R_c synthesis 2- coa_c synthesis independent of pnto_R_c   

two reactions should be synthesized and add to model 
1- MOHMT 3-methyl-2-oxobutanoate hydroxymethyltransferase 3mob_c + h2o_c + mlthf_c → 2dhp_c + thf_c
2-PANTS Pantothenate synthase ala_B_c + atp_c + pant__R_c → amp_c + h_c + pnto__R_c + ppi_c

one reaction shoud be added from reactome 
1- ASP1DC
"""
#%% first reaction pnto_R_c auxotrophy
from cobra import Model, Reaction, Metabolite   

for i in range(len(all_inf)):
    if 'ATCC 334' in all_inf['strain'][i]:
        mod = all_inf['id'][i]+'.json'
        model = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor2/'+mod)
#creating a reaction   
reaction = Reaction('MOHMT')
reaction.name = '3-methyl-2-oxobutanoate hydroxymethyltransferase'
reaction.subsystem = 'Pantothenate and CoA biosynthesis'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default


#adding metabolites
reaction.add_metabolites({
    model.metabolites.get_by_id('3mob_c'): -1,
    model.metabolites.get_by_id('h2o_c'): -1,
    model.metabolites.get_by_id('mlthf_c'): -1,
    model.metabolites.get_by_id('2dhp_c'): 1,
    model.metabolites.get_by_id('thf_c'): 1,
})

#add reaction to model and repair model
reaction.gene_reaction_rule = '(GAP)'
reaction.reaction
model.add_reactions([reaction])
model.repair()

#%% second reaction pnto_R_c auxotrophy
#ala_B_c + atp_c + pant_R_c → amp_c + h_c + pnto_R_c + ppi_c
#creating a reaction   
reaction = Reaction('PANTS')
reaction.name = 'Pantothenate synthase'
reaction.subsystem = 'Pantothenate and CoA biosynthesis'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default


#adding metabolites
reaction.add_metabolites({
    model.metabolites.get_by_id('ala_B_c'): -1,
    model.metabolites.get_by_id('atp_c'): -1,
    model.metabolites.get_by_id('pant_R_c'): -1,
    model.metabolites.get_by_id('amp_c'): 1,
    model.metabolites.get_by_id('h_c'): 1,
    model.metabolites.get_by_id('pnto_R_c'): 1,
    model.metabolites.get_by_id('ppi_c'): 1
})

#add reaction to model and repair model
reaction.gene_reaction_rule = '(GAP)'
reaction.reaction
model.add_reactions([reaction])
model.repair() 

#%% third reaction pnto_R_c auxotrophy


for i in range(len(all_inf)):
    if 'GG (ATCC 53103)' in all_inf['strain'][i]:
        mod = all_inf['id'][i]+'.json'
        template = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor/'+mod)


for i in range(len(all_inf)):
    if 'ATCC 334' in all_inf['strain'][i]:
        mod = all_inf['id'][i]+'.json'
        model = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor2/'+mod)
        
#template = load_json_model('/Users/omidard/Desktop/lacba/ref_model/COBRAModel.json')
target_reactions = []
#rec = ['ASP1DC','PANTS','MOHMT','DPR'] #pnto_R_c for rhamnosus gg
#rec = ['ICDHy','ACONT']#glu solution for gg
#rec = ['IPPMIb','IPPMIa','OMCDC','IPPS','IPMD']#leu solution for gg
rec = ['METS','METB1','HSERTA']# met sol for gg -- template LC2W
#rec = ['ASP1DC'] # pnto_R_c
for i in rec:
    reaction = template.reactions.get_by_id(i)
    reaction2 = Reaction(reaction.id)
    reaction2.name = reaction.name
    reaction2.subsystem = reaction.subsystem
    reaction2.lower_bound = reaction.lower_bound
    reaction2.upper_bound = reaction.upper_bound
    reaction2.add_metabolites(reaction.metabolites)
    reaction2.gene_reaction_rule = '(GAP)'
    target_reactions.append(reaction2)
model.add_reactions(target_reactions)
model.repair()
#%% check if glu solution works
m9(model)
model.reactions.EX_met_L_e.lower_bound = 0
print(model.optimize().fluxes.BIOMASS2)
#%% save model if glu solution worked
cobra.io.json.save_json_model(model,'/Users/omidard/Desktop/lacba/gems/allcor2/'+model.id)
#%% LCW vitamins auxotrophies FP

for i in range(len(all_inf)):
    if 'LC2W' in all_inf['strain'][i]:
        mod = all_inf['id'][i]+'.json'
        model = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor2/'+mod)
        biomass = model.reactions.get_by_id("BIOMASS2")
        #remove metabolites
        #biomass.subtract_metabolites({model.metabolites.get_by_id("thf_c"): -0.00001})
        biomass.subtract_metabolites({model.metabolites.get_by_id("pydx5p_c"): -0.000001})
        #biomass.subtract_metabolites({model.metabolites.get_by_id("btn_c"): -0.00001})
        #biomass.subtract_metabolites({model.metabolites.get_by_id("thmpp_c"): -0.0001})
        #add metabolites
        biomass.add_metabolites({model.metabolites.get_by_id("pydxn_c"): -0.00001})
        biomass.add_metabolites({model.metabolites.get_by_id("ribflv_c"): -0.00001})
        biomass.add_metabolites({model.metabolites.get_by_id("thf_c"): -0.00001})
        biomass.add_metabolites({model.metabolites.get_by_id("btn_c"): -0.00001})
        biomass.add_metabolites({model.metabolites.get_by_id("thm_c"): -0.00001})
        model.repair()
        cobra.io.json.save_json_model(model,'/Users/omidard/Desktop/lacba/gems/allcor2/'+model.id)    

#%% LCW ser_L_c auxotrophies FP step1


sr = []
for i in range(len(all_inf)):
    if 'LC2W' in all_inf['strain'][i]:
        mod = all_inf['id'][i]+'.json'
        model = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor2/'+mod)
        for re in model.genes.GAP.reactions:
            sr.append(re.id)
        for re in model.genes.ORPHAN.reactions:
            sr.append(re.id)
zero = []
for re in sr:
    model = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor2/'+mod)
    m9(model)
    model.reactions.get_by_id(re).lower_bound = 0
    model.reactions.get_by_id(re).upper_bound = 0
    if model.optimize().fluxes.BIOMASS2 == 0:
        zero.append(re)
#%% LCW ser_L_c auxotrophies FP step2
nonz=[]
for i in sr:
    if i not in zero:
        nonz.append(i)
for re in nonz:
    model = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor2/'+mod)
    m9(model)
    model.reactions.get_by_id(re).lower_bound = 0
    model.reactions.get_by_id(re).upper_bound = 0
    model.reactions.EX_ser_L_e.lower_bound = 0
    print(re,' ',model.optimize().fluxes.BIOMASS2 )

#%% second solution for ser FP 
ser = []
for i in range(len(all_inf)):
    if 'LC2W' in all_inf['strain'][i]:
        mod = all_inf['id'][i]+'.json'
        model = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor2/'+mod)
        m9(model)
        flx=model.optimize().fluxes
        for re in model.metabolites.ser_L_c.reactions:
            for i in range(len(flx)):
                if re.id in flx.index[i]:
                    if flx[i] != 0:
                        print(flx.index[i],' ',flx[i], model.reactions.get_by_id(re.id).reaction,' ',model.reactions.get_by_id(re.id).genes)
                        ser.append(re.id)

"""
ser FP solution for LC2W failed
"""
#%%  find possible missing reactions: a general approach         
"""
template = load_json_model('/Users/omidard/Desktop/lacba/ref_model/COBRAModel.json')
m9(template)
template.reactions.EX_leu_L_e.lower_bound = 0
fx1 = template.optimize().fluxes

"""
for i in range(len(all_inf)):
    if 'LC2W' in all_inf['strain'][i]:
        mod = all_inf['id'][i]+'.json'
        template = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor/'+mod)
        m9(template)
        template.reactions.EX_met_L_e.lower_bound = 0
        fx1 = template.optimize().fluxes

for i in range(len(all_inf)):
    if 'GG (ATCC 53103)' in all_inf['strain'][i]:
        mod = all_inf['id'][i]+'.json'
        model = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor2/'+mod)
        m9(model)
        fx2 = model.optimize().fluxes



mis=[]
m2 = [re.id for re in model.reactions]
for re in template.reactions:
    if re.id not in m2:
        if fx1[re.id] != 0:
            mis.append(re.id)
#%% leu solution gg


for i in mis:
    model = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor2/'+mod)
    model.add_reactions([template.reactions.get_by_id(i)])
    m9(model)
    model.reactions.EX_leu_L_e.lower_bound = 0
    if model.optimize().fluxes.BIOMASS2 != 0:
        print(i)
    

# IPPMIb IPPMIa OMCDC IPPS IPMD (leu autothrophy solution for gg -- worked)
#%% methionine autothrophy solution for gg

for i in range(len(all_inf)):
    if 'GG (ATCC 53103)' in all_inf['strain'][i]:
        mod = all_inf['id'][i]+'.json'
        model = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor/'+mod)

gaporf = []
for re in model.genes.GAP.reactions:
    gaporf.append(re.id)
for re in model.genes.ORPHAN.reactions:
    gaporf.append(re.id)
#%%
for re in gaporf:
    model = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor/'+mod)
    m9(model)
    #print('before constraining = ',model.optimize().fluxes.BIOMASS2 )
    model.reactions.get_by_id(re).lower_bound = 0
    model.reactions.get_by_id(re).upper_bound = 0
    if model.optimize().fluxes.BIOMASS2 != 0:
        print(re,' ','thr pos = ', model.optimize().fluxes.BIOMASS2)
        model.reactions.EX_ser_L_e.lower_bound = 0
        print(re,' ','thr neg = ', model.optimize().fluxes.BIOMASS2)
        if model.optimize().fluxes.BIOMASS2 == 0:
            print(re)
    

#%% fine tuning all plantarum growth rates

for i in range(len(all_inf)):
    if 'plantarum' in all_inf['Species'][i]:
        if 'WCFS1' not in all_inf['strain'][i]:
            gaporf=[]
            mod = all_inf['id'][i]+'.json'
            model = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor2/'+mod)
            for re in model.genes.GAP.reactions:
                gaporf.append(re.id)
            for re in model.genes.ORPHAN.reactions:
                gaporf.append(re.id)
            
            for re in gaporf:
                m9(model)
                low = model.reactions.get_by_id(re).lower_bound
                up = model.reactions.get_by_id(re).upper_bound 
                model.reactions.get_by_id(re).lower_bound = 0
                model.reactions.get_by_id(re).upper_bound = 0
                if model.optimize().fluxes.BIOMASS2 > 0.2:
                    print(re,' = ',model.optimize().fluxes.BIOMASS2)
                    reaction = model.reactions.get_by_id(re)
                    reaction.remove_from_model()
                if model.optimize().fluxes.BIOMASS2 <= 0.2:
                    model.reactions.get_by_id(re).lower_bound = low
                    model.reactions.get_by_id(re).upper_bound = up
                cobra.io.json.save_json_model(model,'/Users/omidard/Desktop/lacba/gems/allcor2/'+model.id)
            
#%%
boundc=[]
for i in range(len(all_inf)):
    if 'plantarum' in all_inf['Species'][i]:
        if 'WCFS1' not in all_inf['strain'][i]:
            gaporf=[]
            low = []
            up = []
            mod = all_inf['id'][i]+'.json'
            model = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor/'+mod)
            for re in model.genes.GAP.reactions:
                gaporf.append(re.id)
            for re in model.genes.ORPHAN.reactions:
                gaporf.append(re.id)
            for x in gaporf:
                low.append(model.reactions.get_by_id(x).lower_bound)
                up.append(model.reactions.get_by_id(x).upper_bound)
            bounds = pd.DataFrame()
            bounds['reaction'] = gaporf
            bounds['up'] = up
            bounds['low'] = low
            bounds['model'] = [model.id for j in range(len(gaporf))]
            boundc.append(bounds)
#%%
for inf in boundc[78:610]:
    for j in range(len(inf['reaction'])):
        model = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor2/'+inf['model'][0])
        m9(model)
        model.reactions.get_by_id(inf['reaction'][j]).lower_bound = 0
        model.reactions.get_by_id(inf['reaction'][j]).upper_bound = 0
        
        if model.optimize().fluxes.BIOMASS2 > 0.2:
            print(inf['reaction'][j],' = ',model.optimize().fluxes.BIOMASS2)
            reaction = model.reactions.get_by_id(inf['reaction'][j])
            reaction.remove_from_model()
        if model.optimize().fluxes.BIOMASS2 <= 0.2:
            model.reactions.get_by_id(inf['reaction'][j]).lower_bound = inf['low'][j]
            model.reactions.get_by_id(inf['reaction'][j]).upper_bound = inf['up'][j]
        cobra.io.json.save_json_model(model,'/Users/omidard/Desktop/lacba/gems/allcor2/'+model.id)
        
        
#%% removing unnescessary gaps
for i in range(len(all_inf)):
    if 'plantarum' in all_inf['Species'][i]:
        if 'WCFS1' not in all_inf['strain'][i]:
            gaporf=[]
            mod = all_inf['id'][i]+'.json'
            model = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor/'+mod)
            for re in model.genes.GAP.reactions:
                gaporf.append(re)
            for re in model.genes.ORPHAN.reactions:
                gaporf.append(re)
            
            for r in gaporf:
                r.remove_from_model()
                m9(model)
                print(r.id,' ',model.optimize().fluxes.BIOMASS2)
                if model.optimize().fluxes.BIOMASS2 <= 0.2:
                    model.add_reactions([r])
            cobra.io.json.save_json_model(model,'/Users/omidard/Desktop/lacba/gems/allcor2/'+model.id)
                

#%% growth rate validation

for i in range(len(all_inf)):
    if 'ST-III' in all_inf['strain'][i]:
        mod = all_inf['id'][i]+'.json'
        model = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor/'+mod)
        m9(model)
        gaporf=[]
        print(model.optimize().fluxes.BIOMASS2)
        for re in model.genes.GAP.reactions:
            gaporf.append(re)
        for re in model.genes.ORPHAN.reactions:
            gaporf.append(re)
            
        for r in gaporf:
            r.remove_from_model()
            m9(model)
            print(r.id,' ',model.optimize().fluxes.BIOMASS2)
            if model.optimize().fluxes.BIOMASS2 <= 0.4:
                model.add_reactions([r])
        cobra.io.json.save_json_model(model,'/Users/omidard/Desktop/lacba/gems/allcor2/'+model.id)
            
#%%
for i in range(len(all_inf)):
    if 'DSM 20081' in all_inf['strain'][i]:
        mod = all_inf['id'][i]+'.json'
        model = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor/'+mod)
        m9(model)
        gaporf=[]
        print(model.optimize().fluxes.BIOMASS2)
















