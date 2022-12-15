#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 24 17:46:45 2022

@author: omidard
"""

#%%reactome gapfilling
#1- requirements
import cobra.test
import os
from os.path import join
import pandas as pd
from Bio import SeqIO

#%%2- load reactome model
data_dir = r'/Users/omidard/Desktop'
model = cobra.io.load_matlab_model(join(data_dir, "marlbr.mat"))

#%%3-fba
solution = model.optimize()
print(solution)
model.summary()
#%%
collector =[]
tg = pd.DataFrame(model.genes)
tg.columns = ['id']
tg2 = tg.id[:90000000]
for i in tg2:
    tg3 = str(i)[:10000000000000]
    collector.append(tg3)
tg4 = str(collector)
tg5=tg4.replace("', '","\n")
tg6=tg5.replace(' ','')
tg7=tg6.replace('[','')
tg8=tg7.replace(']','')
tg9=tg8.replace("'","")

reactome = r'/Users/omidard/Desktop/multispecies/Supplementary/prots/reactome.fa'
with open('reactome.fa','w') as of:
    genome = SeqIO.parse(open(reactome),'fasta')
    for seq in genome:
        if seq.id in tg9:
            SeqIO.write([seq], of, "fasta")
