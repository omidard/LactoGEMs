#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  3 18:33:19 2022

@author: omidard


analysis of final models
"""
#imports
import cobra
import pandas as pd
from cobra.io import load_json_model
from glob import glob
import seaborn as sns

#Establish a definition that initializes models to an in silico representation of M9 media

def m9(model):
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
    model.reactions.get_by_id(glc).lower_bound = -30
    model.reactions.get_by_id(arg).lower_bound = -10
    model.reactions.get_by_id(cys).lower_bound = -1
    model.reactions.get_by_id(glu).lower_bound = -10
    model.reactions.get_by_id(ile).lower_bound = -10
    model.reactions.get_by_id(leu).lower_bound = -10
    model.reactions.get_by_id(met).lower_bound = -10
    model.reactions.get_by_id(tyr).lower_bound = -10
    model.reactions.get_by_id(phe).lower_bound = -10
    model.reactions.get_by_id(thr).lower_bound = -10
    model.reactions.get_by_id(val).lower_bound = -10
    model.reactions.get_by_id(cit).lower_bound = -1
    model.reactions.get_by_id(asc).lower_bound = -1
    model.reactions.get_by_id(abz).lower_bound = -10
    model.reactions.get_by_id(ade).lower_bound = -10
    model.reactions.get_by_id(fol).lower_bound = -0.1
    model.reactions.get_by_id(gly).lower_bound = -1
    model.reactions.get_by_id(gua).lower_bound = -10
    model.reactions.get_by_id(ins).lower_bound = -1
    model.reactions.get_by_id(ala).lower_bound = -10
    model.reactions.get_by_id(asp).lower_bound = -2
    model.reactions.get_by_id(his).lower_bound = -10
    model.reactions.get_by_id(lys).lower_bound = -10
    model.reactions.get_by_id(pro).lower_bound = -10
    model.reactions.get_by_id(ser).lower_bound = -1
    model.reactions.get_by_id(trp).lower_bound = -10
    model.reactions.get_by_id(pyd).lower_bound = -0.1
    model.reactions.get_by_id(pdx).lower_bound = -0.1
    model.reactions.get_by_id(rib).lower_bound = -0.1
    model.reactions.get_by_id(ac).lower_bound = -1
    model.reactions.get_by_id(thm).lower_bound = -0.1
    model.reactions.get_by_id(tym).lower_bound = -0.1
    model.reactions.get_by_id(ura).lower_bound = -10
    model.reactions.get_by_id(cl).lower_bound = -1000
    model.reactions.get_by_id(he).lower_bound = -1000
    model.reactions.get_by_id(h2o).lower_bound = -1000
    model.reactions.get_by_id(nh4).lower_bound = -1000
    model.reactions.get_by_id(btn).lower_bound = -0.1
    model.reactions.get_by_id(ca2).lower_bound = -1000
    model.reactions.get_by_id(k).lower_bound = -1000
    model.reactions.get_by_id(co).lower_bound = -1000
    model.reactions.get_by_id(co2).lower_bound = -1000
    model.reactions.get_by_id(pi).lower_bound = -1000
    model.reactions.get_by_id(parametr).upper_bound = 1000
    model.reactions.get_by_id(parametr).lower_bound = 0
    model.reactions.get_by_id(pent).lower_bound = -1
    model.reactions.get_by_id('EX_lac_L_e').upper_bound = 1000
    model.reactions.get_by_id('EX_lac_L_e').lower_bound = 0
    model.objective = 'BIOMASS'
         
   
       
    return model



models2 =glob('%s/*.json'%'Lacticaseibacillusparacasei/gapfilled_models_dir')
allflx=[]
inform = pd.DataFrame()
for ml in models2:
    model = load_json_model(ml)
    m9(model)
    model.reactions.get_by_id('ACAS_2ahbut').lower_bound = -5
    model.reactions.get_by_id('ACAS_2ahbut').upper_bound = 5
    model.reactions.get_by_id('TRSARr').upper_bound = 2
    model.reactions.get_by_id('TRSARr').lower_bound = -2
    model.reactions.get_by_id('P5CRx').upper_bound = 2
    model.reactions.get_by_id('P5CRx').lower_bound = -2
    model.reactions.get_by_id('FRXOXR').upper_bound = 2
    model.reactions.get_by_id('FRXOXR').lower_bound = -2
    model.reactions.get_by_id('DHPR').upper_bound = 2
    model.reactions.get_by_id('DHPR').lower_bound = -2
    model.reactions.get_by_id('ATPM').upper_bound = 5
    model.reactions.get_by_id('ATPM').lower_bound = 5
    model.reactions.get_by_id('PYRt2').upper_bound = 5
    model.reactions.get_by_id('PYRt2').lower_bound = -5
    fba = model.optimize().fluxes
    reactions = list(fba.index)
    selre=[]
    selflx=[]
    for i in reactions:
        if 'EX_' in i:
            selre.append(i)
            flx = fba[i]
            selflx.append(flx)
    inform['Reaction'] = selre
    inform[model.id] = selflx


inf = inform.set_index([pd.Index(selre)])
info = inf.drop(['id'],axis=1)
print(info)



sns.clustermap(info.fillna(0).sort_values(by=info.columns.tolist()), cmap='PiYG', 
              center=0, row_cluster=True, col_cluster=True, figsize=(50,125), vmin=-3, vmax=3)














