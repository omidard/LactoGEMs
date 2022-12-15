#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 11:40:00 2022

@author: omidard
"""

from dgap import all_fluxes, m9,product_associated_genes
import cobra
import pandas as pd
import multiprocessing as mp
from cobra.io import load_json_model
from glob import glob
import numpy as np


models =glob('%s/*.json'%'/home/omidard/gems/allgems/')
pool = mp.Pool(mp.cpu_count())
dfc = pool.map(all_fluxes, [modelid for modelid in models])
allflux = pd.concat(dfc,axis=1)
allflux.fillna(0, inplace=True)
allflux.to_csv('/home/omidard/allflux.csv')
allflux_producers_df = product_associated_genes('EX_mnl_e')
allflux_producers_df.to_csv('/home/omidard/allflux_producers_df.csv')


















    