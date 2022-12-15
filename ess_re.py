#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 13 17:04:36 2022
essential reactions comparison
@author: omidard
"""

from dgap import m9,exchanges_fluxes,total_dir,ess_re
import cobra
import pandas as pd
from cobra.io import load_json_model
from glob import glob
from cobra.manipulation.modify import rename_genes
from cobra import Model, Reaction, Metabolite
import os
import matplotlib.pyplot as plt
import seaborn as sns
import itertools
import math
import multiprocessing as mp


dirs = total_dir('/home/omidard/allgems')
all_genera_essential_genes_list=[]
for dec in dirs:
    directory = dec+'/corr'
    genera = dec.replace('/home/omidard/allgems/','')
    species = dec.replace('/home/omidard/allgems/','')
    models =glob('%s/*.json'%directory)
    species_ess=[]
    for mod in models:
        pool = mp.Pool(mp.cpu_count())
        ess_collector2 = pool.map(ess_re, [mod for mod in models])
        species_ess.append(ess_collector2)
    all_species_essential_genes = pd.concat(species_ess, axis=1)
    all_species_essential_genes.columns = pd.MultiIndex.from_arrays([[genera for i in range(len(all_species_essential_genes.columns))], all_species_essential_genes.columns])
    all_genera_essential_genes_list.append(all_species_essential_genes)
    all_species_essential_genes.to_csv(dec+'/all_species_essential_genes.csv')
all_genera_essential_genes = pd.concat(all_genera_essential_genes_list, axis=1)
all_genera_essential_genes.to_csv('/home/omidard/allgems/all_genera_essential_genes.csv') 
                        
