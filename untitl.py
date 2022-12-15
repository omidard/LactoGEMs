#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 19:34:12 2022

@author: omidard
"""

#%%
#import packages needed
import pandas as pd
from glob import glob
from Bio import Entrez, SeqIO
import cobra
import os
from os.path import join
from os.path import isfile, join
from os import listdir
#%%
directory1 =r'/Users/omidard/Desktop/multispecies/Supplementary/genomes'
files = glob('%s/*.gbk'%directory1)
strain_info = []
for file in files:
        handle = open(file)
        record = list(SeqIO.parse(handle, "genbank"))
        for i in record:
            for f in i.features:
                if f.type=='source':
                    info = {}
                    info['file'] = file.replace('/Users/omidard/Desktop/multispecies/Supplementary/genomes/','')
                    info['id'] = file.replace('/Users/omidard/Desktop/multispecies/Supplementary/genomes/','')
                    for q in f.qualifiers.keys():
                        info[q] = '|'.join(f.qualifiers[q])
                        strain_info.append(info)
                        sinf = pd.DataFrame(strain_info)
                        print(sinf)
#%%
import numpy as np
import biotite
import biotite.sequence as seq
import biotite.sequence.io.fasta as fasta
import biotite.sequence.align as align
import biotite.database.entrez as entrez

# Read each sequence in the file as 'ProteinSequence' object
fasta_file = fasta.FastaFile()
fasta_file.read('targetseqsrxn00011.fasta')
sequences = list(fasta.get_sequences(fasta_file).values())

# BLOSUM62
substitution_matrix = align.SubstitutionMatrix.std_protein_matrix()
# Matrix that will be filled with pairwise sequence identities
identities = np.ones((len(sequences), len(sequences)))
# Iterate over sequences
for i in range(len(sequences)):
    for j in range(i):
        # Align sequences pairwise
        alignment = align.align_optimal(
            sequences[i], sequences[j], substitution_matrix
        )[0]
        # Calculate pairwise sequence identities and fill matrix
        identity = align.get_sequence_identity(alignment)
        identities[i,j] = identity
        identities[j,i] = identity

print(identities)