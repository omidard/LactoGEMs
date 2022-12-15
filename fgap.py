#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  8 04:04:20 2022

@author: omidard
"""

import cobra
import pandas as pd
import seaborn as sns
from cobra.io import load_json_model
from glob import glob
from cobra.manipulation.modify import rename_genes


def m9(model):
    model.solver = 'gurobi'

    for reaction in model.reactions:
        if 'EX_' in  reaction.id:
            reaction.lower_bound=0 
    model.reactions.EX_glc_D_e.lower_bound = -25
    model.reactions.EX_arg_L_e.lower_bound = -10
    model.reactions.EX_cys_L_e.lower_bound = -10
    model.reactions.EX_orot_e.lower_bound = -10
    model.reactions.EX_xan_e.lower_bound = -10
    model.reactions.EX_glu_L_e.lower_bound = -10
    model.reactions.EX_ile_L_e.lower_bound = -10
    model.reactions.EX_leu_L_e.lower_bound = -10
    model.reactions.EX_met_L_e.lower_bound = -10
    model.reactions.EX_tyr_L_e.lower_bound = -10
    model.reactions.EX_phe_L_e.lower_bound = -10
    model.reactions.EX_thr_L_e.lower_bound = -10
    model.reactions.EX_val_L_e.lower_bound = -10
    model.reactions.EX_cit_e.lower_bound = -10
    model.reactions.EX_ascb_L_e.lower_bound = -10
    model.reactions.EX_4abz_e.lower_bound = -10
    model.reactions.EX_ade_e.lower_bound = -10
    model.reactions.EX_fol_e.lower_bound = -0.1
    model.reactions.EX_gly_e.lower_bound = -10
    model.reactions.EX_gua_e.lower_bound = -10
    model.reactions.EX_ins_e.lower_bound = -10
    model.reactions.EX_ala_L_e.lower_bound = -10
    model.reactions.EX_asp_L_e.lower_bound = -10
    model.reactions.EX_his_L_e.lower_bound = -10
    model.reactions.EX_lys_L_e.lower_bound = -10
    model.reactions.EX_pro_L_e.lower_bound = -10
    model.reactions.EX_ser_L_e.lower_bound = -10
    model.reactions.EX_trp_L_e.lower_bound = -10
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
    model.reactions.EX_pnto_R_e.lower_bound = -1
    model.reactions.EX_lac_L_e.upper_bound = 1000
    model.reactions.EX_lac_L_e.lower_bound =0
    model.reactions.ATPM.lower_bound = 1
    model.reactions.ATPM.upper_bound = 1
    model.objective = 'BIOMASS'
    return model


def gapfill_multi(model, missing_genes, **kwargs):
    if 'lower_bound' in kwargs.keys():
        lower_bound = kwargs['lower_bound']
    else:
        lower_bound = model.optimize().objective_value*0.5
    biomass_reactions = [rx.id for rx in model.reactions if rx.objective_coefficient == 1]
    if 'BIOMASS' in kwargs.keys():
        biomass = kwargs['BIOMASS']
        if len(biomass_reactions) > 1:
            for rx in set(biomass_reactions) - {biomass}:
                model.reactions.get_by_id(rx).objective_coefficient = 0        
    else:
        if len(biomass_reactions) > 1:
            raise Exception("This model has more than one objective. \n Please adjust the objective coefficient to 1 for the chosen objective reaction (e.g. biomass or ATP) and 0 for the rest of the reactions, \n or specify the reaction ID to use as an objective.")
        if len(biomass_reactions) > 1:
            raise Exception("The model doesn't have an objective function. Please set the appropriate objective coefficient to 1, or specify the reaction ID to use as an objective.")
        biomass = biomass_reactions[0]
    model.solver.configuration.tolerances.feasibility = 1e-9
    constraints = []
    indicators = []
    for rx in cobra.manipulation.find_gene_knockout_reactions(model, missing_genes):
        indicator = model.problem.Variable('%s_i'%rx.id , type = 'binary')
        indicators.append(indicator)
        new_cstr1 = model.problem.Constraint( rx.flux_expression - rx.upper_bound*indicator ,ub = 0)
        new_cstr2 = model.problem.Constraint(-rx.flux_expression + rx.lower_bound*indicator ,ub = 0)
        constraints += [new_cstr1, new_cstr2]
        model.add_cons_vars([new_cstr1, new_cstr2, indicator])
    model.reactions.get_by_id(biomass).lower_bound = lower_bound
    model.objective = model.problem.Objective(-sum(indicators))
    sol = model.optimize()
    indicator_results = [ind.name[:-2] for ind in indicators if ind.primal != 0.0]
    model.remove_cons_vars(constraints+indicators)
    for rx in set(biomass_reactions):
        model.reactions.get_by_id(rx).objective_coefficient = 1   
    return indicator_results


def missed_genes(orto_matrix,strain_name):
    hom_matrix=pd.read_csv(orto_matrix)
    hom_matrix=hom_matrix.set_index('Unnamed: 0')
    strain=hom_matrix[strain_name]
    print(strain)
    missingGenes=list(strain[strain==0.0].index)
    missingGenes
    return missingGenes


def gaptrack(base_model_dir,orto_matrix,strain_name,target_model_dir,gapfilled_model_dir):
    base=cobra.io.read_sbml_model(base_model_dir)
    m9(base)
    missingGenes = missed_genes(orto_matrix,strain_name)
    gaps = gapfill_multi(base, missingGenes)
    gap=[]
    for re in base.reactions:
        if re.id in gaps:
            gap.append(re)
    print(gap)
    model=cobra.io.load_json_model(target_model_dir)
    model.add_reactions(gap)
    allgen = []
    for i in range(len(gap)):
        for x in gap[i].genes:
            if x.id not in ['EXCHANGE','GAP','BIOMASS','SINK','DEMAND','Diffusion','ORPHAN','spontaneous']:
                allgen.append(x.id)
    rename=[]
    for i in range(len(allgen)):
        rename.append('GAP')  
        rename_dict = {}
    for key in allgen:
        for value in rename:
            rename_dict[key] = value
            rename.remove(value)
            break
    rename_dict
    rename_genes(model,rename_dict)
    cobra.io.save_json_model(model, gapfilled_model_dir+model.id+'.json')
    return model.optimize()



base_model_dir = '/home/omidard/ref_model_dir/marlbr2.xml'
orto_matrix = '/home/omidard/present_absence_dir/ortho_matrix.csv'
target_models = glob('%s/*.json'%'/home/omidard/output_models_dir')
gapfilled_model_dir = '/home/omidard/gapfilled/'
for target_model_dir in target_models:
    strain_nam = target_model_dir.replace('/home/omidard/output_models_dir/','')
    strain_name = strain_nam.replace('.json.json','')
    gaptrack(base_model_dir,orto_matrix,strain_name,target_model_dir,gapfilled_model_dir)





























