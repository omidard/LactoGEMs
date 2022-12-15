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
from cobra.io import load_json_model
from glob import glob
from cobra.manipulation.modify import rename_genes
import defs
import concurrent.futures






gapfill_multi(model, missing_genes)


with concurrent.futures.ThreadPoolExecutor() as executor:
    model_files=glob('%s/*.json'%model_dir)
    models_name = []
    for m in range(len(model_files)):
        mname =model_files[m].replace(model_dir,'')
        modeln = mname.replace('.json','')
        models_name.append(modeln)
    results = executor.map(gapfill_multi, models_name,missing_genes)
    
    



















