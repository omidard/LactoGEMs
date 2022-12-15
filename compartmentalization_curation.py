#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  3 19:32:02 2022

@author: omidard


correction on [c]/[e] compartmentalization from cobratoolbox for compatibility in cobrapy

"""

model=cobra.io.read_sbml_model('')
model.solver = "glpk"

import cobra
import pandas as pd
from cobra.io import load_json_model
from glob import glob
from cobra.manipulation.delete import delete_model_genes, remove_genes
import os

##for removing [] from metabolites compartments

for met in model.metabolites:
    mid = met.id
    if '[c]' in mid or '[e]' in mid:
        mid = mid.replace('[c]','')
        mid = mid.replace('[e]','')
        met.id = mid
        model.repair()
cobra.io.write_sbml_model(model, 'cormet.xml')
