#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 18:16:50 2022

@author: omidard
"""
#%%
#gapfilling based on reactome
import cobra
import os
from os import listdir
from os.path import isfile, join
import subprocess
import pandas as pd
from glob import glob
from Bio import Entrez, SeqIO
from cobra.io import load_json_model

#%%
def get_file_names (directory1,directory2):
    onlyfiles = [f for f in listdir(directory1) if isfile(join(directory1, f))]
    gl = str(onlyfiles)
    gl1=gl.replace("', '","\n")
    gl2=gl1.replace(' ','')
    gl3=gl2.replace('[','')
    gl4=gl3.replace(']','')
    gl5=gl4.replace("'","")
    gl6=gl5.replace(".gbk","")
    gl7=gl6.replace('.DS_Store','')
    with open(directory2+'/targ.txt','w') as tf:
        tf.writelines(gl7)
        tf.close()
    gl8 = []
    txtf = directory2+'/targ.txt'
    with open(txtf) as f:
        for line in f:
            line = line.strip()
            if line != "":
                gl8.append(line)
    gl9=pd.DataFrame({'models':gl8})
    path = directory2+'/modelInformation.xlsx'
    modelInformation = gl9.to_excel(path)
#%%
directory1 =r'/Users/omidard/Desktop/multispecies/Supplementary/Models2'
directory2 =r'/Users/omidard/Desktop/multispecies/Supplementary'
get_file_names (directory1,directory2)
#%%
modelOfInterest=pd.read_excel(directory2+'/modelInformation.xlsx')
modelOfInterest.pop('Unnamed: 0')
print(modelOfInterest)
#%%
def cdm_fba(model,obj):
    for reaction in model.reactions:
        if 'EX_' in  reaction.id:
            reaction.lower_bound=0 
    
    glc = 'EX_glc_D_e'
    arg = 'EX_arg_L_e'
    cys = 'EX_cys_L_e'
    glu = 'EX_glu_L_e'
    ile = 'EX_ile_L_e'
    leu = 'EX_leu_L_e'
    met = 'EX_met_L_e'
    phe = 'EX_phe_L_e'
    thr = 'EX_thr_L_e'
    tyr = 'EX_tyr_L_e'
    val = 'EX_val_L_e'
    cit = 'EX_cit_e'
    asc = 'EX_ascb_L_e'
    abz = 'EX_4abz_e'
    ade = 'EX_ade_e'
    fol = 'EX_fol_e'
    gly = 'EX_gly_e'
    gua = 'EX_gua_e'
    ins = 'EX_ins_e'
    ala = 'EX_ala_L_e'
    asp = 'EX_asp_L_e'
    his = 'EX_his_L_e'
    lys = 'EX_lys_L_e'
    pro = 'EX_pro_L_e'
    ser = 'EX_ser_L_e'
    trp = 'EX_trp_L_e'
    pyd = 'EX_pydam_e'
    pdx = 'EX_pydxn_e'
    rib = 'EX_ribflv_e'
    ac = 'EX_ac_e'
    thm = 'EX_thm_e'
    tym = 'EX_thymd_e'
    ura = 'EX_ura_e'
    cl = 'EX_cl_e'
    he = 'EX_h_e'
    h2o = 'EX_h2o_e'
    nh4 = 'EX_nh4_e'
    btn = 'EX_btn_e'
    ca2 = 'EX_ca2_e'
    k = 'EX_k_e'
    co = 'EX_co_e'
    co2 = 'EX_co2_e'
    pi = 'EX_pi_e'
    model.reactions.get_by_id(glc).lower_bound = -25
    model.reactions.get_by_id(arg).lower_bound = -10
    model.reactions.get_by_id(cys).lower_bound = -10
    model.reactions.get_by_id(glu).lower_bound = -10
    model.reactions.get_by_id(ile).lower_bound = -10
    model.reactions.get_by_id(leu).lower_bound = -10
    model.reactions.get_by_id(met).lower_bound = -10
    model.reactions.get_by_id(tyr).lower_bound = -10
    model.reactions.get_by_id(phe).lower_bound = -10
    model.reactions.get_by_id(thr).lower_bound = -10
    model.reactions.get_by_id(val).lower_bound = -10
    model.reactions.get_by_id(cit).lower_bound = -10
    model.reactions.get_by_id(asc).lower_bound = -10
    model.reactions.get_by_id(abz).lower_bound = -10
    model.reactions.get_by_id(ade).lower_bound = -10
    model.reactions.get_by_id(fol).lower_bound = -10
    model.reactions.get_by_id(gly).lower_bound = -10
    model.reactions.get_by_id(gua).lower_bound = -10
    model.reactions.get_by_id(ins).lower_bound = -10
    model.reactions.get_by_id(ala).lower_bound = -10
    model.reactions.get_by_id(asp).lower_bound = -10
    model.reactions.get_by_id(his).lower_bound = -10
    model.reactions.get_by_id(lys).lower_bound = -10
    model.reactions.get_by_id(pro).lower_bound = -10
    model.reactions.get_by_id(ser).lower_bound = -10
    model.reactions.get_by_id(trp).lower_bound = -10
    model.reactions.get_by_id(pyd).lower_bound = -1000
    model.reactions.get_by_id(pdx).lower_bound = -1000
    model.reactions.get_by_id(rib).lower_bound = -10
    model.reactions.get_by_id(ac).lower_bound = -10
    model.reactions.get_by_id(thm).lower_bound = -1000
    model.reactions.get_by_id(tym).lower_bound = -1000
    model.reactions.get_by_id(ura).lower_bound = -10
    model.reactions.get_by_id(cl).lower_bound = -1000
    model.reactions.get_by_id(he).lower_bound = -1000
    model.reactions.get_by_id(h2o).lower_bound = -1000
    model.reactions.get_by_id(nh4).lower_bound = -1000
    model.reactions.get_by_id(btn).lower_bound = -1000
    model.reactions.get_by_id(ca2).lower_bound = -1000
    model.reactions.get_by_id(k).lower_bound = -1000
    model.reactions.get_by_id(co).lower_bound = -1000
    model.reactions.get_by_id(co2).lower_bound = -1000
    model.reactions.get_by_id(pi).lower_bound = -1000
    model.objective = obj
    solution1 = model.optimize()
    model.summary()
    sol = str(solution1).split(None, 1)[1]
    solution = sol.split(None, 1)[0]
    return print(solution)
#%%list_maker / a list of metabolites wich exist in biomass reaction
def met_ext (model,directory):
    biom = model.reactions.get_by_id('BIOMASS')
    biom2 = str((biom.metabolites).keys())
    bi=biom2.replace('dict_keys(','')
    bi2=bi.replace('Metabolite','\n')
    bi3=bi2.replace('[<','')
    bi4=bi3.replace('[c]','')
    with open(directory+'/biomass_mets.txt','w') as text:
        text.writelines(bi4)
        metlist = []
        metext = directory+'/biomass_mets.txt'
        with open(metext) as f:
            for line in f:
                line = line.strip()
                if line != "":
                    metlist.append(line)
                    metlist3 = []
                    for metab in metlist:
                        metlist2=metab.split(None, 1)[0]
                        metlist3.append(metlist2)
    return pd.DataFrame({'biomass_metabolites':metlist3})
#%%
def rxn_ext(target, model, directory):
    rlist = model.metabolites.get_by_id(target+'[c]').reactions
    rl = str(rlist)
    rl2 = rl.replace('frozenset({','')
    rl3=rl2.replace('>})','')
    rl4 = rl3.replace('<Reaction','\n')
    with open(directory+'/rxnlist.txt','w') as rli:
        rli.writelines(rl4)
        rxnlist = []
        rxntxt = directory+'/rxnlist.txt'
        with open (rxntxt) as f2:
            for line in f2:
                line=line.strip()
                if line != "":
                    rxnlist.append(line)
                    rxnlist3 = []
                    for rxns in rxnlist:
                        rxnlist2 = rxns.split(None, 1)[0]
                        rxnlist3.append(rxnlist2)
    return pd.DataFrame({'target_reactions':rxnlist3})
print (pd.DataFrame({'target_reactions':rxnlist3}))
                    

#%%
def sink_add (model,directory):
    met_ext(model,directory)
    for target in metlist4.biomass_metabolites:
        model.add_boundary(model.metabolites.get_by_id(target+'[c]'), type="sink")
        obj = 'SK_'+target+'[c]'
        cdm_fba(model,obj)
#%% load reactome model
data_dir1 = r'/Users/omidard/Desktop'
rmodel = cobra.io.load_matlab_model(join(data_dir1, "LBR3.mat"))
rmodel_reactions = pd.DataFrame(rmodel.reactions)
#%% iterate through models
def gap_find_fill (tmodel,rmodel):
#%%
for model in modelOfInterest.models:
    missing_reactions = []
    tmodel= cobra.io.load_json_model(directory1+'/'+model)
    tmodel_reactions = pd.DataFrame(tmodel.reactions)
    for reaction in rmodel.reactions:
        if reaction not in tmodel.reactions:
            missing_reactions.append(reaction)
            cdm_fba(tmodel,'BIOMASS')
            if solution == 0:
                met_ext(tmodel,directory)
                for target in metlist4.biomass_metabolites:
                    tmodel2 = tmodel.add_boundary(tmodel.metabolites.get_by_id(target+'[c]'), type="sink")
                    obj = 'SK_'+target+'[c]'
                    cdm_fba(tmodel2,obj)
                    if solution > 0:
                        tmodel3 = tmodel2
                    else:
                        rxn_ext(target, tmodel2, directory)
                       
                        
#%%                  
 modelOfInterest=pd.read_excel(directory2+'/modelInformation.xlsx')
 print(modelOfInterest)                   
                    
#%%
import cobra.test
from cobra.flux_analysis import gapfill
#universal = cobra.io.load_matlab_model(join('/Users/omidard/Desktop/', "LBR3.mat"))

model = load_json_model(r'/Users/omidard/Desktop/multispecies/Supplementary/Models2/GCF_000011985.1.jsonnew.json')
for reaction in model.reactions:
    if 'EX_' in  reaction.id:
        reaction.lower_bound=0 
        
parametr = "BIOMASS"
glc = 'EX_glc_D_e'
arg = 'EX_arg_L_e'
cys = 'EX_cys_L_e'
glu = 'EX_glu_L_e'
ile = 'EX_ile_L_e'
leu = 'EX_leu_L_e'
met = 'EX_met_L_e'
phe = 'EX_phe_L_e'
thr = 'EX_thr_L_e'
tyr = 'EX_tyr_L_e'
val = 'EX_val_L_e'
cit = 'EX_cit_e'
asc = 'EX_ascb_L_e'
abz = 'EX_4abz_e'
ade = 'EX_ade_e'
fol = 'EX_fol_e'
gly = 'EX_gly_e'
gua = 'EX_gua_e'
ins = 'EX_ins_e'
ala = 'EX_ala_L_e'
asp = 'EX_asp_L_e'
his = 'EX_his_L_e'
lys = 'EX_lys_L_e'
pro = 'EX_pro_L_e'
ser = 'EX_ser_L_e'
trp = 'EX_trp_L_e'
pyd = 'EX_pydam_e'
pdx = 'EX_pydxn_e'
rib = 'EX_ribflv_e'
ac = 'EX_ac_e'
thm = 'EX_thm_e'
tym = 'EX_thymd_e'
ura = 'EX_ura_e'
cl = 'EX_cl_e'
he = 'EX_h_e'
h2o = 'EX_h2o_e'
nh4 = 'EX_nh4_e'
btn = 'EX_btn_e'
ca2 = 'EX_ca2_e'
k = 'EX_k_e'
co = 'EX_co_e'
co2 = 'EX_co2_e'
pi = 'EX_pi_e'
model.reactions.get_by_id(glc).lower_bound = -25
model.reactions.get_by_id(arg).lower_bound = -10
model.reactions.get_by_id(cys).lower_bound = -10
model.reactions.get_by_id(glu).lower_bound = -10
model.reactions.get_by_id(ile).lower_bound = -10
model.reactions.get_by_id(leu).lower_bound = -10
model.reactions.get_by_id(met).lower_bound = -10
model.reactions.get_by_id(tyr).lower_bound = -10
model.reactions.get_by_id(phe).lower_bound = -10
model.reactions.get_by_id(thr).lower_bound = -10
model.reactions.get_by_id(val).lower_bound = -10
model.reactions.get_by_id(cit).lower_bound = -10
model.reactions.get_by_id(asc).lower_bound = -10
model.reactions.get_by_id(abz).lower_bound = -10
model.reactions.get_by_id(ade).lower_bound = -10
model.reactions.get_by_id(fol).lower_bound = -10
model.reactions.get_by_id(gly).lower_bound = -10
model.reactions.get_by_id(gua).lower_bound = -10
model.reactions.get_by_id(ins).lower_bound = -10
model.reactions.get_by_id(ala).lower_bound = -10
model.reactions.get_by_id(asp).lower_bound = -10
model.reactions.get_by_id(his).lower_bound = -10
model.reactions.get_by_id(lys).lower_bound = -10
model.reactions.get_by_id(pro).lower_bound = -10
model.reactions.get_by_id(ser).lower_bound = -10
model.reactions.get_by_id(trp).lower_bound = -10
model.reactions.get_by_id(pyd).lower_bound = -10
model.reactions.get_by_id(pdx).lower_bound = -10
model.reactions.get_by_id(rib).lower_bound = -10
model.reactions.get_by_id(ac).lower_bound = -10
model.reactions.get_by_id(thm).lower_bound = -10
model.reactions.get_by_id(tym).lower_bound = -10
model.reactions.get_by_id(ura).lower_bound = -10
model.reactions.get_by_id(cl).lower_bound = -10
model.reactions.get_by_id(he).lower_bound = -10
model.reactions.get_by_id(h2o).lower_bound = -10
model.reactions.get_by_id(nh4).lower_bound = -10
model.reactions.get_by_id(btn).lower_bound = -10
model.reactions.get_by_id(ca2).lower_bound = -10
model.reactions.get_by_id(k).lower_bound = -10
model.reactions.get_by_id(co).lower_bound = -10
model.reactions.get_by_id(co2).lower_bound = -10
model.reactions.get_by_id(pi).lower_bound = -10 
model.reactions.get_by_id(parametr).upper_bound = 10
model.reactions.get_by_id(parametr).lower_bound = -0
model.objective = parametr
solution = model.optimize()
print(solution)
model.summary()

#%% grouping models based on number of reactions
re300_400 =[]
re401_500=[]
re501_600=[]
re601_700=[]
re701_800=[]
re801_900=[]
re901_1000=[]
re1001_1100=[]
re1101_1200=[]
re1201_1300=[]
re1300_1400=[]
re1400_1500 =[]
re1500_1600 =[]
re1600_1700 =[]
models2 =glob('%s/*.json'%'output_models_dir')
for ml in models2:
    ml2 = ml.replace('.json.json','.json')
    model = load_json_model(r'/Users/omidard/Desktop/multispecies/Supplementary/Models2/'+ml2)
    if len(model.reactions) > 300 and len(model.reactions) < 400:
        re300_400.append(model.id)
    if len(model.reactions) > 400 and len(model.reactions) < 500:
        re401_500.append(model.id)
    if len(model.reactions) > 500 and len(model.reactions) < 600:
        re501_600.append(model.id)
    if len(model.reactions) > 600 and len(model.reactions) < 700:
        re601_700.append(model.id)
    if len(model.reactions) > 700 and len(model.reactions) < 800:
        re701_800.append(model.id)
    if len(model.reactions) > 800 and len(model.reactions) < 900:
        re801_900.append(model.id)
    if len(model.reactions) > 900 and len(model.reactions) < 1000:
        re901_1000.append(model.id)
    if len(model.reactions) > 1000 and len(model.reactions) < 1100:
        re1001_1100.append(model.id)
    if len(model.reactions) > 1100 and len(model.reactions) < 1200:
        re1101_1200.append(model.id)
    if len(model.reactions) > 1200 and len(model.reactions) < 1300:
        re1201_1300.append(model.id)
    if len(model.reactions) > 1300 and len(model.reactions) < 1400:
        re1201_1300.append(model.id)
    if len(model.reactions) > 1400 and len(model.reactions) < 1500:
        re1201_1300.append(model.id)
    if len(model.reactions) > 1500 and len(model.reactions) < 1600:
        re1201_1300.append(model.id)
     if len(model.reactions) > 1600 and len(model.reactions) < 1700:
         re1201_1300.append(model.id)
#%%

import pandas as pd
from glob import glob
from Bio import Entrez, SeqIO
import cobra
import os
from os.path import join
from os.path import isfile, join
from os import listdir   
#%%group gapfilling
models=glob('%s/*.json'%'Models2')

parametr = "BIOMASS"
glc = 'EX_glc_D_e'
arg = 'EX_arg_L_e'
cys = 'EX_cys_L_e'
glu = 'EX_glu_L_e'
ile = 'EX_ile_L_e'
leu = 'EX_leu_L_e'
met = 'EX_met_L_e'
phe = 'EX_phe_L_e'
thr = 'EX_thr_L_e'
tyr = 'EX_tyr_L_e'
val = 'EX_val_L_e'
cit = 'EX_cit_e'
asc = 'EX_ascb_L_e'
abz = 'EX_4abz_e'
ade = 'EX_ade_e'
fol = 'EX_fol_e'
gly = 'EX_gly_e'
gua = 'EX_gua_e'
ins = 'EX_ins_e'
ala = 'EX_ala_L_e'
asp = 'EX_asp_L_e'
his = 'EX_his_L_e'
lys = 'EX_lys_L_e'
pro = 'EX_pro_L_e'
ser = 'EX_ser_L_e'
trp = 'EX_trp_L_e'
pyd = 'EX_pydam_e'
pdx = 'EX_pydxn_e'
rib = 'EX_ribflv_e'
ac = 'EX_ac_e'
thm = 'EX_thm_e'
tym = 'EX_thymd_e'
ura = 'EX_ura_e'
cl = 'EX_cl_e'
he = 'EX_h_e'
h2o = 'EX_h2o_e'
nh4 = 'EX_nh4_e'
btn = 'EX_btn_e'
ca2 = 'EX_ca2_e'
k = 'EX_k_e'
co = 'EX_co_e'
co2 = 'EX_co2_e'
pi = 'EX_pi_e'
universal = cobra.io.load_matlab_model(join('/Users/omidard/Desktop/mymat', "marlbr.mat"))
for i in models:
    model = load_json_model(i)
    for reaction in model.reactions:
        if 'EX_' in  reaction.id:
            reaction.lower_bound=0 
    model.solver = 'glpk'
    model.reactions.get_by_id(glc).lower_bound = -30
    model.reactions.get_by_id(arg).lower_bound = -10
    model.reactions.get_by_id(cys).lower_bound = -10
    model.reactions.get_by_id(glu).lower_bound = -10
    model.reactions.get_by_id(ile).lower_bound = -10
    model.reactions.get_by_id(leu).lower_bound = -10
    model.reactions.get_by_id(met).lower_bound = -10
    model.reactions.get_by_id(tyr).lower_bound = -10
    model.reactions.get_by_id(phe).lower_bound = -10
    model.reactions.get_by_id(thr).lower_bound = -10
    model.reactions.get_by_id(val).lower_bound = -10
    model.reactions.get_by_id(cit).lower_bound = -10
    model.reactions.get_by_id(asc).lower_bound = -10
    model.reactions.get_by_id(abz).lower_bound = -10
    model.reactions.get_by_id(ade).lower_bound = -10
    model.reactions.get_by_id(fol).lower_bound = -10
    model.reactions.get_by_id(gly).lower_bound = -10
    model.reactions.get_by_id(gua).lower_bound = -10
    model.reactions.get_by_id(ins).lower_bound = -10
    model.reactions.get_by_id(ala).lower_bound = -10
    model.reactions.get_by_id(asp).lower_bound = -10
    model.reactions.get_by_id(his).lower_bound = -10
    model.reactions.get_by_id(lys).lower_bound = -10
    model.reactions.get_by_id(pro).lower_bound = -10
    model.reactions.get_by_id(ser).lower_bound = -10
    model.reactions.get_by_id(trp).lower_bound = -10
    model.reactions.get_by_id(pyd).lower_bound = -1000
    model.reactions.get_by_id(pdx).lower_bound = -1000
    model.reactions.get_by_id(rib).lower_bound = -10
    model.reactions.get_by_id(ac).lower_bound = -10
    model.reactions.get_by_id(thm).lower_bound = -1000
    model.reactions.get_by_id(tym).lower_bound = -1000
    model.reactions.get_by_id(ura).lower_bound = -10
    model.reactions.get_by_id(cl).lower_bound = -1000
    model.reactions.get_by_id(he).lower_bound = -1000
    model.reactions.get_by_id(h2o).lower_bound = -1000
    model.reactions.get_by_id(nh4).lower_bound = -1000
    model.reactions.get_by_id(btn).lower_bound = -1000
    model.reactions.get_by_id(ca2).lower_bound = -1000
    model.reactions.get_by_id(k).lower_bound = -1000
    model.reactions.get_by_id(co).lower_bound = -1000
    model.reactions.get_by_id(co2).lower_bound = -1000
    model.reactions.get_by_id(pi).lower_bound = -1000
    model.reactions.get_by_id(parametr).upper_bound = 1000
    model.reactions.get_by_id(parametr).lower_bound = -0
    model.objective = parametr
    solution = model.optimize()
    print(solution)
    gapfiller = cobra.flux_analysis.gapfilling.GapFiller(model, universal, demand_reactions=True, exchange_reactions=True, integer_threshold=1e-9)
    gapfiller.model.solver.configuration.tolerances.feasibility = 1e-9
    gapfiller.model.solver.configuration.tolerances.integrality = 1e-9
    miss_re = gapfiller.fill()
    re2 =[]
    for re in miss_re[0]:
        re2.append(re)
        print(len(re2))
    model.add_reactions(re2)
    cobra.io.save_json_model(model, r'/Users/omidard/Desktop/multispecies/Supplementary/Models2/'+model.id+'new.json')
#%%
solution = gapfill(model, universal, demand_reactions=False)
for reaction in solution[0]:
    print(reaction.id)

#%%
import gurobipy
gurobipy.Model()
import optlang.gurobi_interface as grb
grb.Model()
model.solver = 'gurobi'
result = gapfill(model, universal, demand_reactions=False, iterations=20)
for i, entries in enumerate(result):
    print("---- Run %d ----" % (i + 1))
    for e in entries:
        print(e.id)
#%%

import gurobipy
gurobipy.Model()
import optlang.gurobi_interface as grb
grb.Model()
#%%
model.solver = 'glpk'
gapfiller = cobra.flux_analysis.gapfilling.GapFiller(model, universal, demand_reactions=False, integer_threshold=1e-20)
gapfiller.model.solver.configuration.tolerances.feasibility = 1e-20
gapfiller.model.solver.configuration.tolerances.integrality = 1e-20
#gapfiller.model.solver.configuration.tolerances.optimality = 1e-8
miss_re = gapfiller.fill()
#%%
print(len(miss_re[0]))
re2 =[]
for re in miss_re[0]:
    re2.append(re)
model.add_reactions(re2)
#%%













