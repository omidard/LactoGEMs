#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 21:51:55 2022

@author: omidard
"""


import cobra
import pandas as pd
from cobra.io import load_json_model
from glob import glob
from cobra.manipulation.modify import rename_genes
import matplotlib.pyplot as plt
from cobra import Model, Reaction, Metabolite



def m9(model):
    model.solver = 'gurobi'
    for reaction in model.reactions:
        if 'EX_' in  reaction.id:
            reaction.lower_bound=0 
    model.reactions.EX_glc_D_e.lower_bound = -20
    model.reactions.EX_arg_L_e.lower_bound = -1
    model.reactions.EX_cys_L_e.lower_bound = -1
    model.reactions.EX_orot_e.lower_bound = -1
    model.reactions.EX_xan_e.lower_bound = -1
    model.reactions.EX_glu_L_e.lower_bound = -1
    model.reactions.EX_ile_L_e.lower_bound = -1
    model.reactions.EX_leu_L_e.lower_bound = -1
    model.reactions.EX_met_L_e.lower_bound = -1
    model.reactions.EX_tyr_L_e.lower_bound = -1
    model.reactions.EX_phe_L_e.lower_bound = -1
    model.reactions.EX_thr_L_e.lower_bound = -1
    model.reactions.EX_val_L_e.lower_bound = -1
    model.reactions.EX_cit_e.lower_bound = -1
    model.reactions.EX_ascb_L_e.lower_bound = -0.5
    model.reactions.EX_4abz_e.lower_bound = -1
    model.reactions.EX_ade_e.lower_bound = -1
    model.reactions.EX_fol_e.lower_bound = -0.1
    model.reactions.EX_gly_e.lower_bound = -1
    model.reactions.EX_gua_e.lower_bound = -1
    model.reactions.EX_ins_e.lower_bound = -1
    model.reactions.EX_ala_L_e.lower_bound = -1
    model.reactions.EX_asp_L_e.lower_bound = -1
    model.reactions.EX_his_L_e.lower_bound = -1
    model.reactions.EX_lys_L_e.lower_bound = -1
    model.reactions.EX_pro_L_e.lower_bound = -1
    model.reactions.EX_ser_L_e.lower_bound = -1
    model.reactions.EX_trp_L_e.lower_bound = -1
    model.reactions.EX_pydam_e.lower_bound = -0.1
    model.reactions.EX_pydxn_e.lower_bound = -0.1
    model.reactions.EX_ribflv_e.lower_bound = -0.1
    model.reactions.EX_ac_e.lower_bound = -2
    model.reactions.EX_thm_e.lower_bound = -0.1
    model.reactions.EX_thymd_e.lower_bound = -0.1
    model.reactions.EX_ura_e.lower_bound = -2
    model.reactions.EX_cl_e.lower_bound = -1000
    model.reactions.EX_h_e.lower_bound = -1000
    model.reactions.EX_h2o_e.lower_bound = -1000
    model.reactions.EX_nh4_e.lower_bound = -1000
    model.reactions.EX_btn_e.lower_bound = -0.1
    model.reactions.EX_ca2_e.lower_bound = -1000
    model.reactions.EX_k_e.lower_bound = -1000
    model.reactions.EX_co_e.lower_bound = -1000
    model.reactions.EX_co2_e.lower_bound = -1000
    model.reactions.EX_pi_e.lower_bound = -1000
    model.reactions.EX_pnto_R_e.lower_bound = -0.5
    model.reactions.EX_lac_L_e.upper_bound = 1000
    model.reactions.EX_lac_L_e.lower_bound =0
    model.reactions.ATPM.upper_bound = 2
    model.reactions.ATPM.lower_bound = 2
    model.reactions.BIOMASS2.upper_bound = 1000
    model.reactions.BIOMASS2.lower_bound = 0
    model.reactions.BIOMASS.upper_bound = 0
    model.reactions.BIOMASS.lower_bound = 0
    model.objective = 'BIOMASS2'
    return model



dir1='/home/omidard/allgems/inhouse/biomassed2'
dir2='/home/omidard/allgems/inhouse/information2.xls'

models =glob('%s/*.json'%dir1)
gr=[]
mid=[]
reactions_total = []
reactions_gap = []
failed=[] #failed models
temps=[] #gapfilled mopdels
for mod in models:
    model=load_json_model(mod)
    m9(model)
    fba=model.optimize()
    gr.append(fba.fluxes.BIOMASS2)
    mid.append(model.id)
    reactions_gap.append(len(model.genes.GAP.reactions))
    reactions_total.append(len(model.reactions))
all_gems = pd.DataFrame()
all_gems['id']=mid
all_gems['growth']=gr
all_gems['total_reactions']=reactions_total
all_gems['gapfilled_reactions']=reactions_gap

for i in range(len(all_gems.growth)):
    if all_gems.growth[i] < 0.100000:
        failed.append(all_gems.id[i])
    if all_gems.growth[i] > 0.10:
        temps.append(all_gems.id[i])
all_gems.to_excel(dir2)
print(all_gems)
print('failed GEMs = ',len(failed))




