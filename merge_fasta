#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 29 01:38:06 2021

@author: omidard
"""

def omid_ard(input):
    if input % 3 == 0 and input % 5 ==0:
        return "omidard"
    if input % 3 == 0:
        return "omid"
    if input % 5 == 0:
        return "ard"
    return input


print(omid_ard(5))


#%% merge fasta files

import os

DIR = '/Users/omidard/Desktop/fna'

oh = open('one_fasta_file.fasta', 'w')
for f in os.listdir(DIR):
    fh = open(os.path.join(DIR, f))
    for line in fh:
        oh.write(line)
    fh.close()
oh.close()

#%% merge fasta files
#run in command line
cd /pathfile/
cat *.ffn* > newfile.fa

#%%reactions with all possible genes
#load the excell file and grouping genes based on thheir reactions
gpr = pd.read_excel ('gpro.xls') #xls list of gprs and related rxn_ids
gprdf = pd.DataFrame(gpr)
re_ge = dict(tuple(gprdf.groupby('id')))
#concatenate all genes in all strains for each reaction
wholegenes = []
keys = sorted(re_ge.keys())
for each in keys:
    allgenes = list(re_ge.get(each).gpr)
    sent_str = ""
    for i in allgenes:
        sent_str += str(i) + "-"
    sent_str = sent_str[:-1]
    allgenes1 = sent_str.replace('-', 'or')
    allgenes2 = allgenes1.replace(')or(', ') or (')
    allgenes3 = allgenes2.replace('or', ' or ')
    allgenes4 = allgenes3.replace('  or  ', ' or ')
    allgenes5 = allgenes4.replace('(', '')
    allgenes6 = allgenes5.replace(')', '')
    allgenes7 = "({0})".format(allgenes6)
    wholegenes.append(allgenes7)

#make a sorted list out of reactions id
keys = sorted(re_ge.keys())

#save results in a excell file    
#%%reactome plots

import matplotlib.pyplot as plt

plt.bar(range(len(re_ge)), list(re_ge.values()), align='center')
plt.xticks(range(len(re_ge)), list(re_ge.keys()))
#%%
import numpy as np
import matplotlib.pyplot as plt

size = []
for i in range(len(r_g)):
    freq= len(r_g[i])
    size.append(freq)
    
#%% plot
import plotly.graph_objects as go
import plotly.express as px
import numpy as np

x0 = keys
x1 = size

fig = go.Figure()
fig.add_trace(go.Histogram(x=x0))
fig.add_trace(go.Histogram(x=x1))

# The two histograms are drawn on top of another
fig.update_layout(barmode='stack')
fig.show()

#%%
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import figure
plt.style.use('grayscale')

# make data:
np.random.seed(3)
x = keys
y = size

# plot
fig, ax = plt.subplots()

ax.bar(x, y, width=1, edgecolor="white", linewidth=0.5)

ax.set(xlim=(0, 800), xticks=np.arange(1, 800),
       ylim=(0, 50), yticks=np.arange(1, 50))
plt.show()

#%%
import altair as alt
from vega_datasets import data

source = data.wheat()

alt.Chart(source).mark_bar().encode(
    x='year:O',
    y="wheat:Q",
    # The highlight will be set on the result of a conditional statement
    color=alt.condition(
        alt.datum.year == 1810,  # If the year is 1810 this test returns True,
        alt.value('orange'),     # which sets the bar orange.
        alt.value('steelblue')   # And if it's not true it sets the bar steelblue.
    )
).properties(width=600)


#%%
import pandas as pd
gpr = pd.read_excel ('7-jan-LbReactome2_2.xlsx') #xls list of gprs and related rxn_ids
gprdf = pd.DataFrame(gpr)
re_ge = dict(tuple(gprdf.groupby('Subsystem')))                                                            
sorted_rege = sorted(re_ge.values(), key=len)









