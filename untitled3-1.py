#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  6 15:09:46 2021

@author: omidard
"""
#%%find identical seqs
f = open("targetseqsrxn00002.fasta","r")
f1 = open(r'identical.fasta','w')
f2 = open(r'uniq.fasta','w')

seq1=[]
ID=[]

i=0
for line in f:
    i+=1
    if i%2 != 0:
        ID.append(line)
    else:
        seq1.append(line)

for i in range(0,len(seq1)):
    if seq1.count(seq1[i])>=2:
        f1.write(ID[i]+seq1[i])
    else:
        f2.write(ID[i]+seq1[i])

f.close()
f1.close()
f2.close()
#%% multiple pairwase alignment

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

#%%
import subprocess
bashCommand = "mmseqs easy-cluster targetseqsrxn00011.fasta clusterRes tmp --min-seq-id 0.55 -c 0.8 --cov-mode 1"
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
output, error = process.communicate()

#%% grouping genes in one reaction
from Bio import SeqIO
import pandas as pd
import re

fasta_file = r'/Users/omidard/Desktop/reactome/pan.faa'
filenamegenerator = list(re_ge.keys())
joker = list(range(len (re_ge)))
directoryin = r'/Users/omidard/Desktop/reactome/groupedgprsgenes/'
directoryout = r'/Users/omidard/Desktop/reactome/reactionsclusters/'
for i in joker:
    clusname = filenamegenerator[i]
    wanted_file = directoryin+'targetseqs'+str(filenamegenerator[i])+'.fasta'
    tehran = list(SeqIO.parse(wanted_file, "fasta"))
    ker=[]
    for i in range(len(tehran)):
        cph=tehran[i].description
        ker.append(cph)
        ods = str(ker)
        hls = re.findall(r'\[.*?\]', ods)
        osl = str(hls)
        alb = osl[osl.find(' ['):]
        arh = alb.replace("', '", "\n")
        lyn = arh.replace("']", "")
        sob = lyn.replace('",', '\n')
        den = sob.replace("'", "")
        esc = den.replace(" [", "[") 
        nwr = esc.replace(']"','') #final gene names
        nor = ods.replace("['", "")
        prg = nor.replace("']","")
        san = re.sub("[\[].*?[\]]", "", prg)
        pav = san.replace("', '", "")
        rav = pav.replace(" ']", "")
        jav = rav.replace(" ", "\n") #final locoustags name
        for line in jav:
            bag = list(jav.splitlines()) #locustag list
        for line in nwr:
            mos = list(nwr.splitlines()) #gene name list
        earth = {'Locustags':bag, 'Gene_name':mos}
        mway = pd.DataFrame(earth)
        amd = list(tuple(mway.groupby('Gene_name')))
        for i in range(len(amd)):
            gmlt = amd[i][1]['Locustags']
            result_file = directoryout+'cluster'+str(clusname)+'_'+str(i)+'.fasta'
            with open(result_file, "w") as mde:
                for seq in tehran:
                    if seq.id in gmlt:
                        SeqIO.write([seq], mde, "fasta")
  




   
    
    




#%%
database = []
filenamegenerator = list(re_ge.keys())
joker = list(range(len (re_ge)))
for i in joker:
    ver1 = str(list(r_g[i].gpr[:50000])).replace(' ', '\r\n')
    ver2 = ver1.replace(",", "")
    ver3 = ver2.replace("[","")
    ver4 = ver3.replace("]","")
    ver5 = ver4.replace("'","")
    database.append(ver5)
    
    
#%%
fasta_file = r'/Users/omidard/Desktop/reactome/pan.faa' # Input fasta file
for i in joker:
    wanted = database[i]
    result_file = 'ztargetseqs'+str(filenamegenerator[i])+'.fasta'
    wanted = set()
    with open(result_file, "w") as g:
        fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
        for seq in fasta_sequences:
            if seq.id in wanted:
                SeqIO.write([seq], g, "fasta")  






