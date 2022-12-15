#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 13:31:36 2022

@author: omidard
"""

import cobra
import pandas as pd
import seaborn as sns
from cobra.io import load_json_model
from glob import glob
from cobra.manipulation.modify import rename_genes
import matplotlib.pyplot as plt
#%%

model = cobra.io.read_sbml_model('yeastGEM.xml')

#%%

import pandas as pd

exchanges =[]
growth=[]
model = cobra.io.read_sbml_model('yeastGEM.xml')
for re in model.reactions:
    if 'exchange' in re.name:
        exchanges.append(re.id)
        model.solver = 'glpk'
        re.lower_bound=-1000
        re.upper_bound=1000
        model.objective = 'r_4041'
        print(model.optimize().fluxes.r_4041)
        growth.append(model.optimize().fluxes.r_4041)
df = pd.DataFrame()
df['exchanges'] = exchanges
df['growth']= growth


#%%
my_range=list(range(1,len(df.index)+1))
fig, ax = plt.subplots(figsize=(12,4))
plt.vlines(x=my_range, ymin=0, ymax=df['growth'], color='#FF005C', alpha=0.2, linewidth=7)
plt.plot(my_range, df['growth'], "o", markersize=8, color='#FF005C', alpha=0.6)
ax.set_xlabel('exchanges', fontsize=12, fontweight='black', color = '#333F4B')
ax.tick_params(axis='both', which='major', labelsize=10)
ax.set_ylabel('growth rate', fontsize=12, fontweight='black', color = '#333F4B')
# change the style of the axis spines
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_bounds((0, 0.5))
ax.spines['bottom'].set_bounds((1, 0.5))
ax.spines['left'].set_position(('outward', 0.5))
ax.spines['bottom'].set_position(('outward', 0.5))
plt.xticks(my_range, df.models, rotation=90, ha='right')
plt.savefig('arg.png', dpi=300, bbox_inches='tight')
#%%%
from cobra.flux_analysis import flux_variability_analysis
fva = flux_variability_analysis(model, ['r_2056','r_4041'])






