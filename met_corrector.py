#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  3 11:30:03 2022

@author: omidard


model metabolites curation script
"""
#imports 
import cobra
from cobra.io import load_json_model


#transfer numbers from two first characters of metabolites to the end of metabolite name
def met_corrector1(model):
    model=cobra.io.load_json_model(model)
    for metabolite in model.metabolites:
        if metabolite.id[0:2].isdigit():
            if '_c' in metabolite.id:
                print('metabolit name before edit = ',metabolite.id)
                miz = len(metabolite.id)
                metab = metabolite.id[2:miz-2]
                digi = metabolite.id[0:2]
                met = metab+digi+'_c'
                metabolite.id = met
                print('metabolit name after edit = ',met)
                model.repair()
            if '_e' in metabolite.id:
                print('metabolit name before edit = ',metabolite.id)
                miz = len(metabolite.id)
                metab = metabolite.id[2:miz-2]
                digi = metabolite.id[0:2]
                met = metab+digi+'_e'
                print('metabolit name after edit = ',met)
                metabolite.id = met
                model.repair()
    cobra.io.json.save_json_model(model,'corected1'+model.id)
    return model
                
#transfer numbers from first character of metabolites to the end of metabolite name               
def met_corrector2(model):
    model=cobra.io.load_json_model(model)
    for metabolite in model.metabolites:
        if metabolite.id[0:1].isdigit():
            if '_c' in metabolite.id:
                print('metabolit name before edit = ',metabolite.id)
                miz = len(metabolite.id)
                metab = metabolite.id[1:miz-2]
                digi = metabolite.id[0:1]
                met = metab+digi+'_c'
                metabolite.id = met
                print('metabolit name after edit = ',met)
                model.repair()
            if '_e' in metabolite.id:
                print('metabolit name before edit = ',metabolite.id)
                miz = len(metabolite.id)
                metab = metabolite.id[1:miz-2]
                digi = metabolite.id[0:1]
                met = metab+digi+'_e'
                print('metabolit name after edit = ',met)
                metabolite.id = met
                model.repair()
    cobra.io.json.save_json_model(model,'corected2'+model.id)
    return model            
        
        