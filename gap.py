#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 23:37:28 2022

@author: omidard
"""
import cobra
import pandas as pd
from cobra.io import load_json_model
from glob import glob
from cobra.manipulation.delete import delete_model_genes, remove_genes
import os
from os.path import isfile, join
from Bio import Entrez, SeqIO
from os import listdir
from cobra.flux_analysis import gapfill


universal = cobra.io.load_matlab_model(join('ref_model_dir', "marlbr2.mat"))
for reaction in universal.reactions:
    if 'EX_' in  reaction.id:
        reaction.lower_bound=0
parametr = 'BIOMASS'
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
pent = 'EX_pnto_R_e'
universal.reactions.get_by_id(glc).lower_bound = -25
universal.reactions.get_by_id(arg).lower_bound = -10
universal.reactions.get_by_id(cys).lower_bound = -10
universal.reactions.get_by_id(glu).lower_bound = -10
universal.reactions.get_by_id(ile).lower_bound = -10
universal.reactions.get_by_id(leu).lower_bound = -10
universal.reactions.get_by_id(met).lower_bound = -10
universal.reactions.get_by_id(tyr).lower_bound = -10
universal.reactions.get_by_id(phe).lower_bound = -10
universal.reactions.get_by_id(thr).lower_bound = -10
universal.reactions.get_by_id(val).lower_bound = -10
universal.reactions.get_by_id(cit).lower_bound = -10
universal.reactions.get_by_id(asc).lower_bound = -10
universal.reactions.get_by_id(abz).lower_bound = -10
universal.reactions.get_by_id(ade).lower_bound = -10
universal.reactions.get_by_id(fol).lower_bound = -10
universal.reactions.get_by_id(gly).lower_bound = -10
universal.reactions.get_by_id(gua).lower_bound = -10
universal.reactions.get_by_id(ins).lower_bound = -10
universal.reactions.get_by_id(ala).lower_bound = -10
universal.reactions.get_by_id(asp).lower_bound = -10
universal.reactions.get_by_id(his).lower_bound = -10
universal.reactions.get_by_id(lys).lower_bound = -10
universal.reactions.get_by_id(pro).lower_bound = -10
universal.reactions.get_by_id(ser).lower_bound = -10
universal.reactions.get_by_id(trp).lower_bound = -10
universal.reactions.get_by_id(pyd).lower_bound = -1000
universal.reactions.get_by_id(pdx).lower_bound = -1000
universal.reactions.get_by_id(rib).lower_bound = -10
universal.reactions.get_by_id(ac).lower_bound = -10
universal.reactions.get_by_id(thm).lower_bound = -1000
universal.reactions.get_by_id(tym).lower_bound = -1000
universal.reactions.get_by_id(ura).lower_bound = -10
universal.reactions.get_by_id(cl).lower_bound = -1000
universal.reactions.get_by_id(he).lower_bound = -1000
universal.reactions.get_by_id(h2o).lower_bound = -1000
universal.reactions.get_by_id(nh4).lower_bound = -1000
universal.reactions.get_by_id(btn).lower_bound = -1000
universal.reactions.get_by_id(ca2).lower_bound = -1000
universal.reactions.get_by_id(k).lower_bound = -1000
universal.reactions.get_by_id(co).lower_bound = -1000
universal.reactions.get_by_id(co2).lower_bound = -1000
universal.reactions.get_by_id(pi).lower_bound = -1000
universal.reactions.get_by_id(parametr).upper_bound = 10
universal.reactions.get_by_id(parametr).lower_bound = -0
universal.reactions.get_by_id(pent).lower_bound = -1000
universal.objective = 'BIOMASS'   
solution1 = universal.optimize()
print(solution1)
print('-----------universal model has been loaded------------')


model_id =[]
model_size =[]
models2 =glob('%s/*.json'%'output_models_dir')
for ml in models2:
    model1 = load_json_model(ml)
    model_id.append(model1.id)
    model_size.append(len(model1.reactions))
models3 = pd.DataFrame({'model_id':model_id,'model_size':model_size})
models3.sort_values(by=['model_size'], ascending=False, inplace=True)
models4 = list(models3.model_id)
for models in models4:
    model = load_json_model('output_models_dir/'+models+'.json')
    model.solver = 'glpk'
    print(model.id+'/number of models reaction ='+len(model.reactions))
    for reaction in model.reactions:
        if 'EX_' in  reaction.id:
            reaction.lower_bound=0
    parametr = 'BIOMASS'
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
    pent = 'EX_pnto_R_e'
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
    model.reactions.get_by_id(parametr).upper_bound = 10
    model.reactions.get_by_id(parametr).lower_bound = -0
    model.reactions.get_by_id(pent).lower_bound = -1000
    model.objective = 'BIOMASS'  
    print('gapfillin for'+model.id+'has been started')     
    
    solution = gapfill(model, universal, exchange_reactions = True, demand_reactions=True)
    for reaction in solution[0]:
        print(reaction.id)
    re2 =[]
    for re in solution[0]:
        re2.append(re)
        print(len(re2))
    model.add_reactions(re2)
    print('gapfillin for'+model.id+'has been finished')
    cobra.io.save_json_model(model, 'gapfilled_models_dir/'+model.id+'.json')
    
    
    
    
    
    
    
    
    
    