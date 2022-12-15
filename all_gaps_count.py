#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 22:34:19 2022
all gaps analysis
@author: omidard
"""

from dgap import m9,exchanges_fluxes,total_dir
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




dirs = total_dir('/home/omidard/allgems')
all_counts = []
for dec in dirs:
    genera=dec.replace("/home/omidard/allgems/","")
    df=pd.read_csv(dec+'/species_all_gaps.csv')
    allgaps=[]
    for i in range(len(df.index)):
        if i != 0:
            rows = df.loc[df.index[i]]
            for v in range(len(rows)):
                if 'Unnamed: 0' not in rows.index[v]:
                    allgaps.append(rows[v])
    df2=pd.DataFrame({'gaps':allgaps})
    gaps_counter = df2.pivot_table(columns=['gaps'], aggfunc='size')
    gap_counts = pd.DataFrame()
    gap_counts[genera]=gaps_counter
    gap_counts.sort_values(by=[genera],inplace=True,ascending=False)
    all_counts.append(gap_counts)
all_gaps_count = pd.concat(all_counts, axis=1)
all_gaps_count.drop(columns=['inhouse'],inplace=True)
all_gaps_count.to_csv('/home/omidard/allgems/all_gaps_count.csv')

#%%


        