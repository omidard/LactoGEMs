#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 23:22:56 2022

@author: omidard



#########   pull out seq variability from PA matrix   #######3



"""

#%%imports
import pandas as pd
from Bio import SeqIO
from Bio import AlignIO
from Bio import Phylo
import os
from glob import glob
#%%
directory1 =r'/Users/omidard/Desktop/multispecies/Supplementary/bbh'
blast_files = glob('%s/*.csv'%directory1)
genome_id=[]
for i in blast_files:
    gid= i.replace('/Users/omidard/Desktop/multispecies/Supplementary/bbh/reactome_vs_','')
    gid2=gid.replace('_parsed.csv','')
    genome_id.append(gid2)
#%%
allbbh=[]
for i in blast_files:
    bbh=pd.read_csv(i)
    delindx = []
    for i in range(len(bbh)):
        if '->' in bbh.BBH[i]:
            delindx.append(i)
        if bbh.PID[i] <50.0:
            delindx.append(i)
    delindx2 = list(dict.fromkeys(delindx))
    bbh2 = bbh.drop(bbh.index[delindx2])
    bbh3 = pd.DataFrame({'gene':bbh2.gene,'subject':bbh2.subject,'PID':bbh2.PID})
    bbh4 = bbh3.sort_values(by=['PID'], ascending=False)
    allbbh.append(bbh4)
allbbh2 =[]
for i in allbbh:
    i.reset_index(drop=True, inplace=True)
    allbbh2.append(i)
genomes_bbh=pd.DataFrame({'genome_id':genome_id,'bbh':allbbh2}) #final dataframe containing genome id and related bbh results against reactome reference genome

#%%
allgenes =[]
for i in range(len(genomes_bbh)):
     x = genomes_bbh.bbh[i].gene
     for gen in x:
         allgenes.append(gen)
allgenes2 = list(dict.fromkeys(allgenes))
#%%
gend = []
subd = []
pidd = []
for g in allgenes2:
    for i in range(len(genomes_bbh)):
        b = genomes_bbh.bbh[i].gene
        c = genomes_bbh.bbh[i].subject
        p = genomes_bbh.bbh[i].PID
        for z in range(len(b)):
            if g in b[z]:
                gend.append(g)
                subd.append(c[z])
                pidd.append(p[z])
wd = pd.DataFrame({'gene':gend,'subject':subd,'PID':pidd})
wd2 = wd.groupby('gene')
#%%
gooz = wd2.get_group('BAABDPKG_00551')

                
#%%
import cobra
from cobra.io import load_json_model

dire='/Users/omidard/Desktop/initial_models_dir1'
dire2='/Users/omidard/Desktop/lacba/output_models_dir'
f1 = glob('%s/*.json'%dire)
f2 = glob('%s/*.json'%dire2)

b1=[]
for i in f1:
    s = i.replace('/Users/omidard/Desktop/initial_models_dir1/','')
    b1.append(s)

b2 =[]
for i in f2:
    c = i.replace('/Users/omidard/Desktop/lacba/output_models_dir/','')
    c1= c.replace('.json.json','.json')
    b2.append(c1)

g = []

for i in b1:
    if i not in b2:
        op = '/Users/omidard/Desktop/initial_models_dir1/'+i
        model = cobra.io.load_json_model(op)
        cobra.io.save_json_model(model, 'initial_models_dir/'+model.id+'.json')
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
