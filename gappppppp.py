#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 12:12:42 2022

@author: omidard
"""
"""
universal model based on essential reactions in reactome
"""
#%% impoprts
import cobra
import os
from os import listdir
from os.path import isfile, join
import subprocess
import pandas as pd
from glob import glob
from Bio import Entrez, SeqIO
from cobra.io import load_json_model


#%%
universal = cobra.io.load_matlab_model(join('/Users/omidard/Desktop/mymat', "marlbr.mat"))
cobra.flux_analysis.variability.find_essential_genes(model, threshold=0.01, processes=None)