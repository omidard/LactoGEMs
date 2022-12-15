#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 15 19:15:46 2022

@author: omidard
"""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 15 09:27:09 2022

@author: omidard
"""

#diverse template gapfilling: dgap

#imports


import cobra
import pandas as pd
from cobra.io import load_json_model
from glob import glob
from cobra.manipulation.modify import rename_genes
import matplotlib.pyplot as plt
from cobra import Model, Reaction, Metabolite
import multiprocessing as mp
import numpy as np
from joblib import Parallel, delayed, parallel_backend
from joblib import load, dump



#formulation of a chemically defined media
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
    model.reactions.BIOMASS2.upper_bound = 1000
    model.reactions.BIOMASS2.lower_bound = 0
    model.objective = 'BIOMASS2'
    return model



#address to directories
dir1='/Users/omidard/Desktop/lacba/allgems/inhouse/failed'
dir2='/Users/omidard/Desktop/lacba/allgems/inhouse/failed'
dir3='/Users/omidard/Desktop/lacba/allgems/inhouse/feasible/'
dir4='/Users/omidard/Desktop/lacba/allgems/inhouse/failed2/'
    


#scanning models and retrive important information
def scan(directory):
    models =glob('%s/*.json'%directory)
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
    return all_gems,failed,temps
 
        
 
#find best template
def tempfind(all_gems,failed,temps):
    tr=[]
    gapr=[]
    for i in range(len(all_gems)):
        if all_gems.id[i] in temps:
            tr.append(all_gems.total_reactions[i])
            gapr.append(all_gems.gapfilled_reactions[i])
    temps_inf = pd.DataFrame()
    temps_inf['id'] = temps
    temps_inf['total_reactions'] = tr
    temps_inf['gapfilled_reactions']=gapr
    temps_inf.sort_values(by='total_reactions', axis=0, ascending=False, inplace=True, kind='quicksort', na_position='last', ignore_index=True, key=None)
    return temps_inf



#find candidate reactions
def gaps(temps_inf,failed):
    temp_re_id=[] #list1
    template = load_json_model(dir2+'/'+temps_inf.id[len(temps_inf)-1])
    for reaction in template.reactions:
        temp_re_id.append(reaction.id)
    allmissing=[]
    for i in failed:
        missing=pd.DataFrame()
        model = load_json_model(dir2+'/'+i)
        failed_re_id=[] #list two
        for reaction in model.reactions:
            failed_re_id.append(reaction.id)
        missing_reactions = list(set(temp_re_id).difference(failed_re_id))
        missing[model.id] = missing_reactions
        allmissing.append(missing)
    return allmissing


#find gaps
def fgap(allmissing):
    all_gaps=[]
    for i in range(len(allmissing)):
        gap=[]
        gap_per_gem = pd.DataFrame()
        for x in allmissing[i][allmissing[i].columns[0]]:
            template = load_json_model(dir2+'/'+temps_inf.id[len(temps_inf)-1])
            m9(template)
            template.reactions.get_by_id(x).lower_bound = 0
            template.reactions.get_by_id(x).upper_bound = 0
            if template.optimize().fluxes.BIOMASS2 == 0:
                gap.append(x)
        gap_per_gem[allmissing[i].columns[0]]=gap
        print(gap_per_gem)
        all_gaps.append(gap_per_gem)
    return all_gaps

#add candidate reactions and save models based on gapfilling status
def addgaps(all_gaps,temps_inf):
    for i in range(len(all_gaps)):
        model = load_json_model(dir2+'/'+all_gaps[i].columns[0])
        template = load_json_model(dir2+'/'+temps_inf.id[len(temps_inf)-1])
        for x in all_gaps[i][all_gaps[i].columns[0]]:
            reaction = template.reactions.get_by_id(x)
            reaction2 = Reaction(reaction.id)
            reaction2.name = reaction.name
            reaction2.subsystem = reaction.subsystem
            reaction2.lower_bound = reaction.lower_bound
            reaction2.upper_bound = reaction.upper_bound
            reaction2.add_metabolites(reaction.metabolites)
            reaction2.gene_reaction_rule = '(GAP)'
            model.add_reactions([reaction2])
            model.repair()
            print(reaction2.id)
            print('-----done')
            m9(model)
            print(model.optimize().fluxes.BIOMASS2)
        if model.optimize().fluxes.BIOMASS2>0.01:
            cobra.io.json.save_json_model(model,dir3+model.id)
        if model.optimize().fluxes.BIOMASS2<=0.01:
            cobra.io.json.save_json_model(model,dir4+model.id)








output = scan(dir1)
all_gems=output[0]
failed=output[1]
temps=output[2]
temps_inf = tempfind(all_gems,failed,temps)
allmissing = gaps(temps_inf,failed)
all_gaps = gapfind(allmissing)







"""
#find gaps
def fgap(allmissing):
    all_gaps=[]
    for i in range(len(allmissing)):
        gap=[]
        gap_per_gem = pd.DataFrame()
        for x in allmissing[i][allmissing[i].columns[0]]:
            template = load_json_model(dir2+'/'+temps_inf.id[len(temps_inf)-1])
            m9(template)
            template.reactions.get_by_id(x).lower_bound = 0
            template.reactions.get_by_id(x).upper_bound = 0
            if template.optimize().fluxes.BIOMASS2 == 0:
                gap.append(x)
        gap_per_gem[allmissing[i].columns[0]]=gap
        print(gap_per_gem)
        all_gaps.append(gap_per_gem)
    return all_gaps

"""

