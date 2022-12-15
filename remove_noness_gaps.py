#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  2 09:49:59 2022

@author: omidard
"""

from dgap import m9,gap_extract,none_ess_gap
import cobra
import pandas as pd
from cobra.io import load_json_model
from glob import glob
from cobra.manipulation.modify import rename_genes
from cobra import Model, Reaction, Metabolite
import os




names = []
rootdir = '/home/omidard/allgems'
for file in os.listdir(rootdir):
    local = os.path.join(rootdir, file)
    if os.path.isdir(local):
        names.append(local)

for name in names:
    os.system("cd "+name)
    os.mkdir('corr')
    directory1 = name+'/feasible/'
    directory2 = name+'/corr/'
    models =glob('%s/*.json'%directory1)
    for mod in models:
        model = load_json_model(mod)
        total_gaps = gap_extract(model)
        none_ess_gap(total_gaps,directory1,directory2)
        
        
    


#%%
import os
import pexpect
names = []
rootdir = '/Users/omidard/Desktop/lacba/allgems/allgems'
for file in os.listdir(rootdir):
    local = os.path.join(rootdir, file)
    server = local.replace('/Users/omidard/Desktop/lacba/allgems','omidard@20.123.39.83:/home/omidard')
    file = server.replace('omidard@20.123.39.83:/home/omidard/allgems/','')
    if os.path.isdir(local):
        names.append(file)
        print('scp -r '+local+'/feasible '+server+'/feasible/')
        """
        child = pexpect.spawn('scp -r '+local+'/feasible '+server+'/feasible/')
        child.expect ('Password:',timeout=None)
        child.sendline ('Biosustain522')
        """
#%%%

import os
names = []
rootdir = '/home/omidard/allgems'
for file in os.listdir(rootdir):
    local = os.path.join(rootdir, file)
    if os.path.isdir(local):
        names.append(local)
        
        
        
        
        
        
        
        
        
        
        
        
        