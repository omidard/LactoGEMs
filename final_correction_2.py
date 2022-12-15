#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 17:25:59 2022

@author: omidard
"""
from cobra.io import load_json_model,save_json_model
from dgap import m9
import cobra
import pandas as pd
from cobra.manipulation.modify import rename_genes
from cobra import Model, Reaction, Metabolite
all_inf=pd.read_csv('/Users/omidard/Desktop/lacba/inf.csv')
#%%
gr=[]
mid=[]
models=[]
from dgap import m9
for mo in all_inf['id']:
    mod = mo+'.json'
    model = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor/'+mod)
    models.append(model)
    m9(model)
    flx = model.optimize().fluxes
    gr.append(flx.BIOMASS2)
    mid.append(model.id)
df = pd.DataFrame()
df['id']=mid
df['growth_rate']=gr

#%%
gr=[]
mid=[]
for model in models:

    #model = models[5]
    m9(model)
    flx = model.optimize().fluxes
    print(model.id) 
    print('-----------------') 
    print('Products')
    for i in range(len(flx)):
        if flx[i] >0:
            if 'EX_' in flx.index[i]:
                print(str(flx.index[i])+' = '+str(flx[i]))
    print('-----------------')            
    print('Substrates')
    for i in range(len(flx)):
        if flx[i] <-0.5:
            if 'EX_' in flx.index[i]:
                print(str(flx.index[i])+' = '+str(flx[i]))
    print('-----------------')
    print('Growth rate')
    print(flx.BIOMASS2)
    print('-----------------')
    
    print('-----------------')
    print('Positive unrealistic fluxes')
    for i in range(len(flx)):
        if flx[i] >10:
            if 'EX_' not in flx.index[i]:
                print(str(flx.index[i])+' = '+str(flx[i]))
    print('-----------------')
    print('Negative unrealistic fluxes')
    for i in range(len(flx)):
        if flx[i] <-10:
            if 'EX_' not in flx.index[i]:
                print(str(flx.index[i])+' = '+str(flx[i]))
    print('#########################################################') 
     
    gr.append(flx.BIOMASS2)
    mid.append(model.id)
    df = pd.DataFrame()
    
df['id']=mid
df['growth_rate']=gr


#%%
weak=[]
zer=[]
for i in range(len(df)):
    if df['growth_rate'][i]<0.1:
        weak.append(df['id'][i])
    if df['growth_rate'][i]==0:
        zer.append(df['id'][i])


#%% zero all reactions
allzerore=[]
for i in zer:
    model = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor/'+i)
    for re in model.reactions:
        if re.id not in allzerore:
            allzerore.append(re.id)

#%%weak all genes reactions
allweakre=[]
for i in weak:
    model = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor/'+i)
    for re in model.reactions:
        if re.id not in allweakre:
            allweakre.append(re.id)
    
#%%missing reactions from zeros
missing_reactions = list(set(allweakre).difference(allzerore))

#%% all weak and zeros
templre = []
templ = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor/GCF_900518845.1.json')
m9(templ)
print(templ.optimize().fluxes.BIOMASS2)
for re in templ.reactions:
    templre.append(re.id)


target = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor/GCF_000712375.1.json')
m9(target)
print(target.optimize().fluxes.BIOMASS2)
targre=[]
for re in target.reactions:
    targre.append(re.id)


misre = list(set(templre).difference(targre))


#%%

for i in templre:
    if templ.reactions.get_by_id(i).lower_bound != target.reactions.get_by_id(i).lower_bound:
        print(i,' = ',' lowerbounds ', templ.reactions.get_by_id(i).lower_bound,target.reactions.get_by_id(i).lower_bound)
        target.reactions.get_by_id(i).lower_bound = templ.reactions.get_by_id(i).lower_bound
    if templ.reactions.get_by_id(i).upper_bound != target.reactions.get_by_id(i).upper_bound:
        print(i,' = ',' upperbounds ', templ.reactions.get_by_id(i).upper_bound,target.reactions.get_by_id(i).upper_bound)
        target.reactions.get_by_id(i).upper_bound = templ.reactions.get_by_id(i).upper_bound



#%% missing reactions which could carry fluxes
reactions_list=[]
m9(templ)
fx = templ.optimize().fluxes
nzfx = []
for i in range(len(fx)):
    if fx[i] != 0:
        nzfx.append(fx.index[i])

m9(target)
print(target.optimize().fluxes.BIOMASS2)
tre=[]
for re in target.reactions:
    tre.append(re.id)
ttt=[]
for i in nzfx:
    if i not in tre:
        ttt.append(i)
        
#%%



for i in range(len(fx.index)):
    for j in misre:
        if j == fx.index[i]:
            if fx[i] != 0:
                print(fx.index[i],' = ',model.reactions.get_by_id(j).reaction,' ',fx[i])
                reactions_list.append(fx.index[i])

#%%
mreactions=[]
ref = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor/GCF_900518845.1.json')
for i in misre:
    reaction = ref.reactions.get_by_id(i)
    reaction2 = Reaction(reaction.id)
    reaction2.name = reaction.name
    reaction2.subsystem = reaction.subsystem
    reaction2.lower_bound = reaction.lower_bound
    reaction2.upper_bound = reaction.upper_bound
    reaction2.add_metabolites(reaction.metabolites)
    reaction2.gene_reaction_rule = '(ORPHAN)'
    mreactions.append(reaction2)
#%% all to all models
#for mo in all_inf['id']:
    #mod = mo+'.json'

target.add_reactions(mreactions)
target.repair()
m9(target)
print(target.optimize().fluxes.BIOMASS2)
    #save_json_model(model,'/Users/omidard/Desktop/lacba/gems/allcor/'+model.id)
#%% add to weak models

for mo in all_inf['id']:
    mod = mo+'.json'
    model = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor/'+mod)
    model.add_reactions(mreactions)
    model.repair()
    save_json_model(model,'/Users/omidard/Desktop/lacba/gems/allcor/'+model.id)

#%% add to zero models
for mod in zer:
    model = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor/'+mod)
    model.add_reactions(mreactions)
    model.repair()
    m9(model)
    print(model.id,' = ',model.optimize().fluxes.BIOMASS2)
    save_json_model(model,'/Users/omidard/Desktop/lacba/gems/allcor/'+model.id)




#%% tuning upper and lower bounds of reactions (no nead to re-run)

for i in range(len(ptar)):
    if ptar['Growth'][i] <= 0.099:
        mod = ptar['id'][i]
        ind = temps_df[temps_df['Species'] == ptar['Species'][i]].index.tolist()
        template_mod = temps_df['id'][ind[0]]
        model = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor/'+mod)
        m9(model)
        template = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor/'+template_mod)
        m9(template)
        # find missing reactions from target model
        temp_re=[]
        for re in template.reactions:
            temp_re.append(re.id)
        targ_re=[]
        for re in model.reactions:
            targ_re.append(re.id)
        
        
        for i in targ_re:
            if i in temp_re:
                if model.reactions.get_by_id(i).lower_bound != template.reactions.get_by_id(i).lower_bound:
                    print(i)
                    model.reactions.get_by_id(i).lower_bound = template.reactions.get_by_id(i).lower_bound
                if model.reactions.get_by_id(i).upper_bound != template.reactions.get_by_id(i).upper_bound:
                    print(i)
                    model.reactions.get_by_id(i).upper_bound = template.reactions.get_by_id(i).upper_bound
        save_json_model(model,'/Users/omidard/Desktop/lacba/gems/allcor2/'+model.id)
#%%
from glob import glob

models =glob('%s/*.json'%'/Users/omidard/Desktop/lacba/gems/allcor2')
for mod in models:
    model=load_json_model(mod)
    m9(model)
    print(model.id,' = ',model.optimize().fluxes.BIOMASS2)






#%%automated detection of best template 
#1- zero model is from which species?
#2- nonzero model from same species as template
spe=[]
for i in range(len(all_inf['id'])):
    if all_inf['id'][i]+'.json' in zer:
        spe.append(all_inf['Species'][i])
sp = []
for i in spe:
    if i not in sp:
        sp.append(i)
        
#%%extract all gems within target species
possible_temps=[]
spp=[]
for i in range(len(all_inf)):
    if all_inf['Species'][i] in sp:
        possible_temps.append(all_inf['id'][i]+'.json')
        spp.append(all_inf['Species'][i])
ptar= pd.DataFrame()
ptar['id']=possible_temps
ptar['Species']=spp

gro=[]
for i in possible_temps:
    model = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor/'+i)
    m9(model)
    gro.append(model.optimize().fluxes.BIOMASS2)
ptar['Growth'] = gro

#%%select highest growth rate per species

ptarg = ptar.groupby(['Species'],sort=True)
spgl = []
for i in ptarg:
    spg = i
    spgl.append(spg)
spgl2 = spgl

sspgl=[]
for i in spgl2:
    s = i[1].sort_values(by=['Growth'],ascending=False,ignore_index=True)
    sspgl.append(s) #list of possible template per target species
#%% 
temps_list = []
for i in sspgl:
    t = i[0:1]
    temps_list.append(t)
#%% templates dataframe
temps_df= pd.concat(temps_list,axis=0) #call templates from this df
temps_df.reset_index(inplace=True)
#%%

#find zero models and their best template and load them
for i in range(len(ptar)):
    if ptar['Growth'][i] <= 0.099:
        mod = ptar['id'][i]
        ind = temps_df[temps_df['Species'] == ptar['Species'][i]].index.tolist()
        template_mod = temps_df['id'][ind[0]]
        model = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor/'+mod)
        m9(model)
        template = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor/'+template_mod)
        m9(template)
        # find missing reactions from target model
        temp_re=[]
        for re in template.reactions:
            temp_re.append(re.id)
        targ_re=[]
        for re in model.reactions:
            targ_re.append(re.id)
        misre = list(set(temp_re).difference(targ_re)) 
        print('missing reactions has been found')
        reactions_list=[]
        fx = template.optimize().fluxes
        for i in range(len(fx.index)):
            for j in misre:
                if j == fx.index[i]:
                    if fx[i] != 0:
                        print(fx.index[i],' = ',template.reactions.get_by_id(j).reaction,' ',fx[i])
                        reactions_list.append(fx.index[i])

        mreactions=[]
        ref = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor/'+template_mod)
        for i in reactions_list:
            reaction = ref.reactions.get_by_id(i)
            reaction2 = Reaction(reaction.id)
            reaction2.name = reaction.name
            reaction2.subsystem = reaction.subsystem
            reaction2.lower_bound = reaction.lower_bound
            reaction2.upper_bound = reaction.upper_bound
            reaction2.add_metabolites(reaction.metabolites)
            reaction2.gene_reaction_rule = '(GAP)'
            mreactions.append(reaction2)
        
        
        model.add_reactions(mreactions)
        model.repair()
        m9(model)
        print('#################>>>>',model.optimize().fluxes.BIOMASS2)
        save_json_model(model,'/Users/omidard/Desktop/lacba/gems/allcor2/'+model.id)

#%%
#%%
#%%
#%%
newzer=[]
for i in zer:
    model = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor/'+i)
    m9(model)
    if model.optimize().fluxes.BIOMASS2 == 0:
        newzer.append(i)
#%%
spe=[]
for i in range(len(all_inf['id'])):
    if all_inf['id'][i]+'.json' in weak:
        spe.append(all_inf['Species'][i])
sp = []
for i in spe:
    if i not in sp:
        sp.append(i)
        
#%%extract all gems within target species
possible_temps=[]
spp=[]
gr=[]
for i in range(len(all_inf)):
    if all_inf['Species'][i] in sp:
        possible_temps.append(all_inf['id'][i]+'.json')
        spp.append(all_inf['Species'][i])
ptar= pd.DataFrame()
ptar['id']=possible_temps
ptar['Species']=spp

for i in ptar['id']:
    ind = df[df['id'] == i].index.tolist()
    gr.append(df['growth_rate'][ind[0]])
ptar['Growth'] = gr

#%%select highest growth rate per species

ptarg = ptar.groupby(['Species'],sort=True)
spgl = []
for i in ptarg:
    spg = i
    spgl.append(spg)
spgl2 = spgl

sspgl=[]
for i in spgl2:
    s = i[1].sort_values(by=['Growth'],ascending=False,ignore_index=True)
    sspgl.append(s) #list of possible template per target species
#%% 
temps_list = []
for i in sspgl:
    t = i[0:1]
    temps_list.append(t)
#%% templates dataframe
temps_df= pd.concat(temps_list,axis=0) #call templates from this df
temps_df.reset_index(inplace=True)

#%% manual template selection
template = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor/GCF_000008065.1.json')
m9(template)
temp_re=[]
for re in template.reactions:
    temp_re.append(re.id)
fx = template.optimize().fluxes


for i in sspgl[1]['id']:
    model = load_json_model('/Users/omidard/Desktop/lacba/gems/allcor/'+i)
    m9(model)
    targ_re=[]
    for re in model.reactions:
        targ_re.append(re.id)
    misre = list(set(temp_re).difference(targ_re)) 
    print('missing reactions has been found')
    reactions_list=[]
    for i in range(len(fx.index)):
        for j in misre:
            if j == fx.index[i]:
                if fx[i] != 0:
                    print(fx.index[i],' = ',template.reactions.get_by_id(j).reaction,' ',fx[i])
                    reactions_list.append(fx.index[i])

    mreactions=[]
    ref = template
    for i in reactions_list:
        reaction = ref.reactions.get_by_id(i)
        reaction2 = Reaction(reaction.id)
        reaction2.name = reaction.name
        reaction2.subsystem = reaction.subsystem
        reaction2.lower_bound = reaction.lower_bound
        reaction2.upper_bound = reaction.upper_bound
        reaction2.add_metabolites(reaction.metabolites)
        reaction2.gene_reaction_rule = '(GAP)'
        mreactions.append(reaction2)
    
    
    model.add_reactions(mreactions)
    model.repair()
    m9(model)
    print('#################>>>>',model.optimize().fluxes.BIOMASS2)
    save_json_model(model,'/Users/omidard/Desktop/lacba/gems/allcor2/'+model.id)

 
            
            
            
            
            
            
            
            
            
            
            
            