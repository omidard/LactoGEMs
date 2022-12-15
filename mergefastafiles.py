#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 19:25:55 2022

@author: omidard
"""

#merge fasta files


cd /pathfile/
cat *.ffn* > newfile.fa


#%%
import os
newfile = 'mergedref.fa'
pathfile = r'/Users/omidard/Desktop/reflactobacillaceae/ncbi_dataset/ncbi_dataset/data/faa'

cmd_line='/Users/omidard/Desktop/reflactobacillaceae/ncbi_dataset/ncbi_dataset/data/faa/cat *.ffn* > newfile.fa'

print ('running merge fasts with following command line...')
print (cmd_line)
os.system(cmd_line)
#%%
from Bio import SeqIO
ref = []
genome = SeqIO.parse(r'/Users/omidard/Desktop/reflactobacillaceae/ncbi_dataset/ncbi_dataset/data/faa/mergedref.fa','fasta')
for seq in genome:
    if seq.id not in ref:
        ref.append(seq.id)
#%%
ref2 = list(dict.fromkeys(ref))
#%%
with open('ref.faa','w') as r:
    genome = SeqIO.parse(r'/Users/omidard/Desktop/reflactobacillaceae/ncbi_dataset/ncbi_dataset/data/faa/mergedref.fa','fasta')
    for seq in genome:
        if seq.id in ref:
            SeqIO.write([seq],r,'fasta')
