#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 16:50:43 2022

@author: omidard
"""

#reduce genome

import pandas as pd
from Bio import SeqIO
from Bio import AlignIO
from Bio import Phylo
import os
#%% extract locus tags from reactome and remove repeated locus tags
reactome = pd.read_excel('mar_lbr_apr.xlsx')
allgprstr = reactome.GPR
gprdf = str(list(allgprstr))
loc1=gprdf.replace("', '","\n")
loc2=loc1.replace("['", "")
loc3=loc2.replace("']", "")
loc4=loc3.replace(" or ", "\n")
loc5 = loc4.replace('(','')
loc6 = loc5.replace(')','')
loc7 = loc6.replace(' and ','\n')
loc8 = list(loc7.splitlines())
loc9 = list(dict.fromkeys(loc8))
#%% pull out those genes frome reactome reference genome and writ it again to ilimitane non metabolic genes. this will help blast to run 10 times faster
genome_file = 'reactome.fa'
result_file = 'rnreactome.fa'
with open(result_file, "w") as g:
    fasta_sequences = SeqIO.parse(open(genome_file),'fasta')
    for seq in fasta_sequences:
        if seq.id in loc9:
            SeqIO.write([seq], g, "fasta")
#%%













