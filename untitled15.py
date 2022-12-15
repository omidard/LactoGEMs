#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 22:02:53 2022

@author: omidard
"""
#find gaps based on biomass precursor connectivities




from dgap import scan,tempfind,gaps,zx,fluxanalyze,addgaps,m9
import multiprocessing as mp
import cobra
import pandas as pd
from cobra.io import load_json_model
from glob import glob
from cobra.manipulation.modify import rename_genes
import matplotlib.pyplot as plt
from cobra import Model, Reaction, Metabolite
import multiprocessing as mp
import numpy as np

model = load_json_model('GCF_902374445.1.json')
temp = load_json_model('allgems/Lactobacillus_iners/biomassed/GCF_009857205.1.json')


def precursor(model):
    prec=[]
    info =[]
    for met in model.reactions.BIOMASS2.metabolites:
        model.add_boundary(model.metabolites.get_by_id(met.id), type="sink")
        m9(model)
        model.objective = 'SK_'+met.id
        re = 'SK_'+met.id
        info.append(met.id)
        info.append(str(model.optimize()))
        fba = model.optimize().fluxes
        for i in range(len(fba)):
            if re in fba.index[i]:
                if fba[i] == 0:
                    prec.append(met.id)
    return prec,info
            
        


def pretargets(template,prec):
    alltarg=[]
    for pre in prec:
        targmet = temp.metabolites.get_by_id(pre)
        targre = []
        for re in targmet.reactions:
            targre.append(re.id)
        alltarg.append(targre)
    return alltarg
        
        
        
        
        
def modelre(model,met):
    targmet = model.metabolites.get_by_id(met)
    modre = []
    for re in targmet.reactions:
        modre.append(re.id)
    return modre





def missing (list1,modre):
    misre = list(set(list1).difference(modre))
    return misre
    
        
        
        
        
        
        
        
        