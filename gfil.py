#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  3 17:39:05 2022

@author: omidard
"""

#gfill
import cobra
import pandas as pd
from cobra.io import load_json_model
from glob import glob
from cobra.manipulation.delete import delete_model_genes, remove_genes
import cobra
import os
from os.path import join
import pandas as pd
from glob import glob
from Bio import Entrez, SeqIO
import cobra
import os
from os.path import join
from os.path import isfile, join
from os import listdir
import os
import pandas as pd
from glob import glob
from Bio import Entrez, SeqIO
import cobra
from os.path import join
from os.path import isfile, join
from os import listdir


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



universal = cobra.io.load_matlab_model(join('/Users/omidard/Desktop/mymat', "marlbr2.mat"))
models2 =glob('%s/*.json'%'/Users/omidard/Desktop/multispecies/Supplementary/Models2')
for lb_models in models2:
    model = cobra.io.load_json_model(lb_models)
    model.solver = 'glpk'
    tr =[]
    for reactions in universal.reactions:
        if reactions not in model.reactions:
            tr.append(reactions)
    modelCopy=model.copy()
    modelCopy.add_reactions(tr)
    er = cobra.flux_analysis.find_essential_reactions(modelCopy, threshold=None, processes=None)
    ftr = []
    for reactions in modelCopy.reactions:
        if reactions in tr and reactions not in er:
            ftr.append(reactions)
            model.sdd_reactions(ftr)
            cobra.io.save_json_model(model, 'gapfilled_models_dir/'+model.id+'gapfill.json')
            
            
            
            






















