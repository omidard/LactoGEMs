#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  3 21:21:46 2022

@author: omidard
"""
import cobra
from cobra.io import load_json_model
from glob import glob
from cobra.manipulation.modify import rename_genes
from cobra import Model, Reaction, Metabolite
from dgap import biomass2
import cobra
import pandas as pd
import multiprocessing as mp
from cobra.io import load_json_model
from glob import glob
import numpy as np

models =glob('%s/*.json'%'/home/omidard/corrected/')
pool = mp.Pool(mp.cpu_count())
pool.map(biomass2, [mod for mod in models])