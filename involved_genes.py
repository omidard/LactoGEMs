#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 10:21:41 2022
product_essential_genes
@author: omidard
"""

from dgap import product_based_essential_genes,m9
import cobra
import pandas as pd
import multiprocessing as mp
from cobra.io import load_json_model
from glob import glob
import numpy as np

models =glob('%s/*.json'%'/home/omidard/gems/allgems/')
pool = mp.Pool(mp.cpu_count())
dfc = pool.map(product_based_essential_genes, [mod for mod in models[0:500]])
inv_genes = pd.concat(dfc,axis=1)
inv_genes.fillna(0, inplace=True)
inv_genes.to_csv('/home/omidard/inv_genes_4abut_e.csv')
    


            
    
    
    

        
                        
        