#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 01:15:38 2022

@author: omidard
"""

from dgap import scan,tempfind,gaps,zx,fluxanalyze,addgaps,m9
import multiprocessing as mp
import cobra
import pandas as pd
from cobra.io import load_json_model
from glob import glob
from dgap import m9

models = glob('%s/*.json'%'/home/omidard/allgems/inhouse/feasible/')
for mod in models:
    model = load_json_model(mod)
    m9(model)
    print(model.optimize().fluxes.BIOMASS2)
    