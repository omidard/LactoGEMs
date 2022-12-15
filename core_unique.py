#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 10:03:54 2022

@author: omidard
"""


from dgap import m9, all_reactions, gene_associated_reactions
import cobra
import pandas as pd
import multiprocessing as mp
from cobra.io import load_json_model
from glob import glob
import numpy as np


models =glob('%s/*.json'%'/home/omidard/gems/allgems/')
pool = mp.Pool(mp.cpu_count())
dfc = pool.map(gene_associated_reactions, [modelid for modelid in models])
core_shell = pd.concat(dfc,axis=1)
core_shell.to_csv('/home/omidard/all_reactions.csv')












