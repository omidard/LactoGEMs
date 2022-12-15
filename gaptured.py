#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  2 15:21:21 2022

@author: omidard

gapture gapfilling
"""

#imports
import cobra
import pandas as pd
import seaborn as sns
from cobra.io import load_json_model
from glob import glob
from cobra.manipulation.modify import rename_genes







#Establish a definition that initializes models to an in silico representation of M9 media

def m9(model):
    model.solver = 'gurobi'
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
    model.reactions.get_by_id(fol).lower_bound = -0.1
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
    model.reactions.get_by_id(pyd).lower_bound = -0.1
    model.reactions.get_by_id(pdx).lower_bound = -0.1
    model.reactions.get_by_id(rib).lower_bound = -0.1
    model.reactions.get_by_id(ac).lower_bound = -10
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
    model.reactions.get_by_id(pent).lower_bound = -1
    model.reactions.get_by_id('BIOMASS').lower_bound = 0
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
        print(biomass)
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
    
    
    # removing changes to model
    model.remove_cons_vars(constraints+indicators)
    for rx in set(biomass_reactions):
        model.reactions.get_by_id(rx).objective_coefficient = 1 
    
        
    return indicator_results



#%%
dir_one = '/Users/omidard/Desktop/lacba/Lacticaseibacillusparacasei/ref_model_dir/marlbr2.xml'
dir_two = '/Users/omidard/Desktop/lacba/Lactobacillus_crispatus/output_models_dir'
dir_three = '/Users/omidard/Desktop/lacba/Lactobacillus_crispatus/present_absence_dir/ortho_matrix.csv'
dir_four = '/Users/omidard/Desktop/lacba/Lactobacillus_crispatus/output_models_dir/'
dir_five = '/Users/omidard/Desktop/lacba/Lactobacillus_crispatus/gapfilled/'

#%%

model_files=glob('%s/*.json'%dir_two)

for m in model_files:
    mname = m.replace(dir_four,'')
    modeln = mname.replace('.json.json','')
    #Gather the list of base strain genes that have no homolog in strain of interest, an input to the below function
    hom_matrix=pd.read_csv(dir_three)
    hom_matrix=hom_matrix.set_index('Unnamed: 0')
    strain=hom_matrix[modeln]
    missingGenes=list(strain[strain==0.0].index)



    model=cobra.io.load_json_model(dir_four+mname)
    m9(model)
    model.optimize()

    
    base=cobra.io.read_sbml_model(dir_one)
    m9(base)
    base.optimize()


    gaps = gapfill_multi(base, missingGenes)
    print('><><><>','gapfilling is finished')
    #print('gaps = ', gaps)
    gaptured = []
    for re in base.reactions:
        if re.id in gaps:
            gaptured.append(re)

        

    

    model.add_reactions(gaptured)
    
    
    allgen = []
    for i in range(len(gaptured)):
        for x in gaptured[i].genes:
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
    #print(rename_dict) 

    rename_genes(model,rename_dict)
    model.repair()
    print('><><><>','missing genes has been renamed to gap')
    m9(model)
    print(model.optimize())
    cobra.io.save_json_model(model, dir_five+model.id+'.json')
    print('><><><>','model reactions gapfilled = ', len(model.reactions))
    print('done')



























