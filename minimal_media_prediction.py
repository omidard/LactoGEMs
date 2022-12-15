#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 11 09:42:57 2022
minimal media exploration
@author: omidard
"""

from dgap import m9,exchange_reactions,essential_substrates, total_dir
import cobra
import pandas as pd
from cobra.io import load_json_model
from glob import glob
from cobra.manipulation.modify import rename_genes
from cobra import Model, Reaction, Metabolite
import os
import multiprocessing as mp

directory = '/home/omidard/gems/allgems/'
substrates = exchange_reactions(directory)
models =glob('%s/*.json'%directory)
substrates = essential_substrates(models,substrates)
substrates[0].to_csv('/home/omidard/minimal_media.csv')
