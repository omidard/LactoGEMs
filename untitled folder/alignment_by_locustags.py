#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 30 14:30:39 2021

@author: omidard
"""

#%% step1: required imports:
import pandas as pd
from Bio import SeqIO
from Bio import AlignIO
from Bio import Phylo
import os


#%% step2: required inputs:
gpr = pd.read_excel ('gpro.xls') #xls list of gprs and related rxn_ids
fasta_file = "pan.faa" # Input fasta file
#%%step3: grouping all genes based on the reactions involved in
gprdf = pd.DataFrame(gpr)
re_ge = dict(tuple(gprdf.groupby('id')))
r_g = list(re_ge.values())
#%%step4: save genes list related to each reaction into a text file (recomended)
filenamegenerator = list(re_ge.keys())
joker = list(range(len (re_ge)))
for i in joker:
    with open ('genes'+str(filenamegenerator[i])+'.txt', "w") as m:
        ver1 = str(list(r_g[i].gpr[:50000])).replace(' ', '\r\n')
        ver2 = ver1.replace(",", "")
        ver3 = ver2.replace("[","")
        ver4 = ver3.replace("]","")
        ver5 = ver4.replace("'","")
        m.writelines(ver5)
        m.close()
#%%step5: retrive fasta seqs for each reaction
fasta_file = "pan.faa" # Input fasta file
for i in joker:
    wanted_file = 'genes'+str(filenamegenerator[i])+'.txt'
    result_file = 'targetseqs'+str(filenamegenerator[i])+'.fasta'
    wanted = set()
    with open(wanted_file) as f:
        for line in f:
            line = line.strip()
            if line != "":
                wanted.add(line)
                with open(result_file, "w") as g:
                    fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
                    for seq in fasta_sequences:
                        if seq.id in wanted:
                            SeqIO.write([seq], g, "fasta")
                            
                            
                        
#%%step6: run MSA
#conda create -n clustalw2 -c biobuilds -y clustalw
clustalw_exe = "/Users/omidard/opt/anaconda3/envs/clustalw2/bin/clustalw2"
for i in joker:
    from Bio.Align.Applications import ClustalwCommandline
    cline = ClustalwCommandline("clustalw2", infile = 'targetseqs'+str(filenamegenerator[i])+'.fasta')
    #print(cline)
    clustalw_cline = ClustalwCommandline(clustalw_exe, infile = 'targetseqs'+str(filenamegenerator[i])+'.fasta')
    assert os.path.isfile(clustalw_exe), "Clustal_W executable is missing"
    stdout, stderr = clustalw_cline()
    #print(clustalw_cline)
    #ClustalAlign = AlignIO.read('targetseqs'+str(filenamegenerator[i])+'.aln', "clustal")                                                      
#%%print alingment result
#from Bio import AlignIO    
ClustalAlign = AlignIO.read("targetseqsrxn00012.aln", "clustal")
print(ClustalAlign)    
#%%prind dendrogram
#from Bio import Phylo
tree = Phylo.read("targetseqsrxn00001.dnd", "newick")
Phylo.draw_ascii(tree)