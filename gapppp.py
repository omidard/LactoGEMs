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
    model.objective = obj
    solution1 = model.optimize()
    print(solution1)
    #model.summary()
    #sol = str(solution1).split(None, 1)[1]
    #solution = sol.split(None, 1)[0]
    #return print(solution)
#%%set parameters
universal = cobra.io.load_matlab_model(join('ref_model_dir', "marlbr2.mat"))
models2 =glob('%s/*.json'%'output_models_dir')
for lb_models in models2:
    model = cobra.io.load_json_model(lb_models)
    print(len(model.reactions))
    for reaction in model.reactions:
        if 'EX_' in  reaction.id:
            reaction.lower_bound=0
            
    
            
            cdm_fba(model,obj)
            gapfiller = cobra.flux_analysis.gapfilling.GapFiller(model, universal, demand_reactions=True, exchange_reactions=True, integer_threshold=1e-200)
            gapfiller.model.solver.configuration.tolerances.feasibility = 1e-100
            gapfiller.model.solver.configuration.tolerances.integrality = 1e-100
            miss_re = gapfiller.fill()
            re2 =[]
            for re in miss_re[0]:
                re2.append(re)
                print(len(re2))
            model.add_reactions(re2)
            cobra.io.save_json_model(model, 'gapfilled_models_dir/'+model.id+'.json')


















