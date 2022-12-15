#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 10:25:44 2022
genome-gems analysis
@author: omidard
"""

from dgap import m9,exchanges_fluxes,total_dir
from glob import glob
import pandas as pd
from Bio import SeqIO
import os
import re
import cobra
import pandas as pd
from cobra.io import load_json_model
import seaborn as sns



#%% genomic information

dirs = total_dir('/Users/omidard/Desktop/lacba/gems/genomes')
alldf = []
for dec in dirs:
    directory1=dec+'/genbank/'
    files = glob('%s/*.gbk'%directory1)
    files2 = sorted(files)
    strain_info = []
    for file in files2:
            handle = open(file)
            record = list(SeqIO.parse(handle, "genbank"))
            trap=[]
            for i in record:
                for f in i.features:
                    if f.type=='source':
                        info = {}
                        info['file'] = file.replace(directory1,'')
                        info['id'] = file.replace(directory1,'')
                        info['genome_size'] = len(i.seq)/1000000
                        info['contigs']=len(record)
                        for q in f.qualifiers.keys():
                            info[q] = '|'.join(f.qualifiers[q])
                            strain_info.append(info)
                    if f.type == 'CDS':
                        if 'hypothetical' in f.qualifiers['product'][0] or 'Hypothetical' in f.qualifiers['product'][0]:
                            trap.append(1)
                    info['hypo']=len(trap)
    sinf = pd.DataFrame(strain_info)
    for i in range(len(sinf.index)):
        sinf['file'][i]=sinf['id'][i].replace('.gbk','')
    sinf.set_index('file',inplace=True)
    sinf2 = sinf.drop_duplicates(subset='id', keep='first', inplace=False, ignore_index=False)
    sinf3 = sinf2.sort_values(by=['id'])
    for i in range(len(sinf3.index)):
        sinf3['id'][i]=sinf3['id'][i].replace('.gbk','')
    alldf.append(sinf3)
    
all_genus_genomes_inf = pd.concat(alldf, axis=0)


#%%curation
directory1 ='/Users/omidard/Desktop/lacba/gems/allgems'
files = glob('%s/*.json'%directory1)
gems=[]
for i in files:
    b=i.replace('/Users/omidard/Desktop/lacba/gems/allgems/','')
    e=b.replace('.json','')
    gems.append(e)
misgems=[]
for i in range(len(all_genus_genomes_inf)):
    if all_genus_genomes_inf['id'][i] not in gems:
        misgems.append(all_genus_genomes_inf['id'][i])


adf = all_genus_genomes_inf.drop(misgems, axis=0)
adf.sort_values(by='id',inplace=True)
aadf = adf.drop_duplicates()




#%%% gems information
import cobra
import pandas as pd
from cobra.io import load_json_model

from dgap import m9,exchanges_fluxes,total_dir

directory2 = '/Users/omidard/Desktop/lacba/gems/allgems'
gaps=[]
genes=[]
reactions=[]
for i in aadf.index:
    mod = directory2+'/'+i+'.json'
    model = load_json_model(mod)
    gaps.append(len(model.genes.GAP.reactions))
    genes.append(len(model.genes))
    reactions.append(len(model.reactions))
aadf['gaps']=gaps
aadf['genes']=genes
aadf['reactions']=reactions

#%%metadata

genomes_inf=pd.read_csv('/Users/omidard/Desktop/lacba/gems/df_antismash.csv')
genomes_inf.set_index('genome_id', inplace=True)
genomes_inf.drop(['source', 'species','strain','closest_placement_reference',
                  'organism','assembly','tax_id','#Seq','Avg','Min','Max',
                  'refseq_category','refseq','genbank','refseq_genbank_identity',
                  'biosample','assembly_type','release_type','genome_representation',
                  'pubMLST','Kingdom','Phylum','Class','Order',
                  'Organism','Genus_cat','bgcs_count','bgcs_on_contig_edge',
                  'protoclusters_count','cand_clusters_count','Total bp','contigs'], 
                 axis=1,inplace=True)
genomes_inf.drop(['GCF_009556455.1'], axis=0)


rowsc=[]
for i in range(len(aadf.index)):
    for v in genomes_inf.index:
        if aadf.index[i] == v:
            rows = genomes_inf.loc[v]
            rowsc.append(rows)
row = pd.concat(rowsc, axis=1)
row2=row.T
all_inf= pd.concat([row2,aadf], axis=1)


#%%plot
sns.lmplot(
    data=all_inf, x="gc_content", y="genome_len", col="Family", hue="Family",
    col_wrap=2, palette="muted", ci=None,
    height=10, scatter_kws={"s": 9, "alpha": 0.75}
)
#%%plot
sns.set_theme(style="darkgrid")
g = sns.jointplot(x="reactions", y="gaps", data=all_inf,
                  kind="scatter", truncate=False,
                  color="red", height=12)

#%%%cool plot
sns.pairplot(all_inf, hue="assembly_level",height=1.5,kind="kde",diag_kind='kde')

#%% heatmap plot                  
                            
sns.set_theme(style="white")

g = sns.JointGrid(data=all_inf, x="reactions", y="genome_len", space=0)
g.plot_joint(sns.kdeplot,
             fill=True,
             thresh=0, levels=100, cmap="rocket")
g.plot_marginals(sns.histplot, color="#03051A", alpha=1, bins=25)
    


#%%% Genomes information plot

import seaborn as sns
sns.set_theme(style="whitegrid")

# Load the example planets dataset

cmap = sns.cubehelix_palette(rot=-.2, as_cmap=True)
g = sns.relplot(
    data=all_inf,
    x="gaps", y="hypo",
    hue="assembly_level", size="Genus", sizes=(10, 400),height=12
)
g.set(xscale="log", yscale="log")
g.ax.xaxis.grid(True, "minor", linewidth=.5)
g.ax.yaxis.grid(True, "minor", linewidth=.5)
g.despine(left=True, bottom=True)


#%%

import seaborn as sns
import matplotlib.pyplot as plt

#sns.set_theme(style="ticks")
sns.set_context("paper", font_scale=1.1) 
# Initialize the figure with a logarithmic x axis
f, ax = plt.subplots(figsize=(4, 8))
#ax.set_xscale("log")


# Plot the orbital period with horizontal boxes
sns.boxplot(x="gc_content", y="Species", data=all_inf.sort_values(by='gc_content'),
            whis=[1, 10], width=0.7, color='#1c6474',saturation=1,fliersize=2,linewidth=0.7)



# Tweak the visual presentation
ax.xaxis.grid(False)
ax.set(ylabel="")
sns.despine(trim=False, left=False,right=False,top=False)


#%%%
from matplotlib import cbook
from matplotlib import cm
from matplotlib.colors import LightSource
import matplotlib.pyplot as plt
import numpy as np

# Load and format data

z = all_inf['contigs']
x = all_inf['auN']
y = all_inf['hypo']


fig = plt.figure(figsize=(15, 8))
ax = fig.add_subplot(projection='3d')


ax.scatter(x, y, z)

ax.set_xlabel('auN')
ax.set_ylabel('Hypothetical proteins')
ax.set_zlabel('contigs')

plt.show()

#%%% number of gems for each species-genera
from collections import Counter
genus = list(all_inf['Genus'])
species = list(all_inf['Species'])
Count_genus = Counter(genus)
Count_species = Counter(species)

f, ax = plt.subplots(figsize=(6, 8))
sns.barplot(x=[Count_species[i] for i in Count_species], y=[i for i in Count_species],
            label="Total", color="#5e8874")
ax.set(xlim=(0, 700), ylabel="",
       xlabel="Number of GEMs per Genus")
sns.despine(left=True, bottom=True)

#%%

import matplotlib.pyplot as plt

fig, ax = plt.subplots()

fruits = [i for i in Count_genus]
counts = [Count_genus[i] for i in Count_genus]

ax.bar(fruits, counts, color="#5e8874")
#%%%
ax.set_ylabel('fruit supply')
ax.set_title('Fruit supply by kind and color')
ax.legend(title='Fruit color')

plt.show()




























