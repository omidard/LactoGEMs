#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  3 14:25:59 2022

@author: omidard
"""

import cobra
from cobra.io import load_json_model




#formulation of a chemically defined media
def cdm(model):
    model.solver = 'gurobi'
    for reaction in model.reactions:
        if 'EX_' in  reaction.id:
            reaction.lower_bound=0

    #amino acids
    model.reactions.EX_arg_L_e.lower_bound = -2
    model.reactions.EX_cys_L_e.lower_bound = -2
    model.reactions.EX_glu_L_e.lower_bound = -2
    model.reactions.EX_ile_L_e.lower_bound = -2
    model.reactions.EX_leu_L_e.lower_bound = -2
    model.reactions.EX_met_L_e.lower_bound = -2
    model.reactions.EX_tyr_L_e.lower_bound = -2
    model.reactions.EX_phe_L_e.lower_bound = -2
    model.reactions.EX_thr_L_e.lower_bound = -2
    model.reactions.EX_val_L_e.lower_bound = -2
    model.reactions.EX_gly_e.lower_bound = -2
    model.reactions.EX_ala_L_e.lower_bound = -2
    model.reactions.EX_asp_L_e.lower_bound = -2
    model.reactions.EX_his_L_e.lower_bound = -2
    model.reactions.EX_lys_L_e.lower_bound = -2
    model.reactions.EX_pro_L_e.lower_bound = -2
    model.reactions.EX_ser_L_e.lower_bound = -2
    model.reactions.EX_trp_L_e.lower_bound = -2


    #carbon
    model.reactions.EX_glc_D_e.lower_bound = -15
    #model.reactions.EX_glc_D_e.upper_bound = -5
    model.reactions.EX_ac_e.lower_bound = -1
    model.reactions.EX_cit_e.lower_bound = -1

    #nuclotides
    model.reactions.EX_thymd_e.lower_bound = -0.1
    model.reactions.EX_ura_e.lower_bound = -1
    model.reactions.EX_gua_e.lower_bound = -1
    model.reactions.EX_ins_e.lower_bound = -1
    model.reactions.EX_ade_e.lower_bound = -1
    model.reactions.EX_xan_e.lower_bound = -1
    model.reactions.EX_orot_e.lower_bound = -1
    
    #vitamins
    model.reactions.EX_btn_e.lower_bound = -0.1
    model.reactions.EX_pnto_R_e.lower_bound = -0.5
    model.reactions.EX_thm_e.lower_bound = -0.1
    model.reactions.EX_pydam_e.lower_bound = -0.1
    model.reactions.EX_pydxn_e.lower_bound = -0.1
    model.reactions.EX_ribflv_e.lower_bound = -0.1
    model.reactions.EX_fol_e.lower_bound = -0.1
    model.reactions.EX_ascb_L_e.lower_bound = -0.5
    model.reactions.EX_4abz_e.lower_bound = -1
    model.reactions.EX_nac_e.lower_bound = -1
    
    #minerals
    #model.reactions.EX_o2_e.lower_bound = -10
    model.reactions.EX_cl_e.lower_bound = -10
    model.reactions.EX_h_e.lower_bound = -1000
    model.reactions.EX_h2o_e.lower_bound = -10
    model.reactions.EX_nh4_e.lower_bound = -10
    model.reactions.EX_ca2_e.lower_bound = -10
    model.reactions.EX_co_e.lower_bound = -10
    model.reactions.EX_co2_e.lower_bound = -10
    model.reactions.EX_pi_e.lower_bound = -100
    model.reactions.EX_cobalt2_e.lower_bound = -10
    model.reactions.EX_cu2_e.lower_bound = -10
    model.reactions.EX_fe3_e.lower_bound = -10
    model.reactions.EX_k_e.lower_bound = -10
    model.reactions.EX_mn2_e.lower_bound = -10
    model.reactions.EX_so4_e.lower_bound = -10
    model.reactions.EX_na1_e.lower_bound = -10
    model.reactions.EX_mg2_e.lower_bound = -10
    model.reactions.EX_zn2_e.lower_bound = -10
   
    #correction bounds HCLTr
    for reaction in model.reactions:
        if 'HCLTr' in  reaction.id:
            reaction.lower_bound=-10
            reaction.upper_bound=10
        if 'PTRCt2' in  reaction.id:
            reaction.lower_bound=-10
    #model.reactions.MNLpts.lower_bound = -1
    model.reactions.ATPM.upper_bound = 1
    model.reactions.ATPM.lower_bound = 1 #0.38
    model.reactions.ATPM.upper_bound = 1 #0.38
    model.reactions.EX_ptrc_e.upper_bound = 0
    model.reactions.EX_ptrc_e.lower_bound = 0
    #model.reactions.PTRCORNt7.bounds=(-10,10)
    
    
    #model.reactions.BIOMASS2.upper_bound = 1000
    #model.reactions.BIOMASS2.lower_bound = 0
    model.objective = 'BIOMASS2'
    return model