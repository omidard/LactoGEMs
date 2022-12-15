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
from Bio.Align.Applications import ClustalwCommandline
#%%
def read_value_from_excel(filename, column="C", row=1):

    return pd.read_excel(filename, skiprows=row - 1, usecols=column, nrows=1, header=None, names=["Value"]).iloc[0]["Value"]
#%% run the function to get the gene list for each reaction
row = 318 #change the numerical value based on reaction's row in excell

fasta_file = "/Users/omidard/Desktop/multispecies/ref/reactome.fa" # Input fasta file
react = r'/Users/omidard/Desktop/7-jan-LbReactome2_2.xlsx' #excell file of reactome
reaction_genes = read_value_from_excel (react, 'C',row)
print('\n','\n','-----------------------------------------------------------------','\n','locus tags','\n','-----------------------------------------------------------------','\n')
#% find gene ids for manual boolean rule curation
#reaction_genes = 'AADJCKDE_00878','ADJCKDE_01855','AADJCKDE_01854'
gene_description = []
involved_genes =[]
input_genome2 = r'/Users/omidard/Desktop/OTHER/faa/reactome.faa' #merged faa file of 49 lactobacilli genome used for reactome reconstruction
genome2 = SeqIO.parse(open(input_genome2),'fasta')
for seq in genome2:
    if seq.id in reaction_genes:
        gene_description.append(seq.description)

locustags = [s.split(None, 1)[0] for s in gene_description]
short_list = [s.split(None, 1)[1] for s in gene_description]
short_list2 = pd.DataFrame(short_list)
short_list2.columns = ['enzymes']
short_list3 = list(short_list2.groupby('enzymes'))
for i in range(len(short_list3)):
    desc = short_list3[i][0]
    involved_genes.append(desc)
statement = str(involved_genes)
statement2 = statement.replace("['",'\n')
statement3 = statement2.replace("', '","\n")
statement4 = statement3.replace("']","\n")
print (locustags)
print ('\n','\n','------number of sequences =',len(locustags),'---------')
print('\n','#--------------------------------------------------------------#')
print(' #--------------------reaction information----------------------#')
print(' #--------------------------------------------------------------#')
print('\n','>>>>>>>>>','Modelseed_id =',read_value_from_excel (react, 'O', row))
print('\n','>>>>>>>>>','Reaction =',read_value_from_excel (react, 'B', row))
print('\n','>>>>>>>>>',"Reaction's name = ",read_value_from_excel (react, 'D', row))
print('\n','-----------------------------------------------------------------','\n','>>>>> different type of genes has been found in this reaction are <<<<<<<','\n','-----------------------------------------------------------------','\n',statement4,'\n','-----------------------------------------------------------------','\n')
#%%
fasta_sequence = SeqIO.parse(open(fasta_file),'fasta')
sequences = r'/Users/omidard/Desktop/alignment/seqs.fasta'
with open(sequences, "w") as g:
    for seq in fasta_sequence:
        if seq.id in locustags:
            SeqIO.write([seq], g, "fasta")
#%%step6: run MSA
#conda create -n clustalw2 -c biobuilds -y clustalw
clustalw_exe = "/Users/omidard/opt/anaconda3/bin/clustalw2"
cline = ClustalwCommandline("clustalw2", infile = r'/Users/omidard/Desktop/alignment/seqs.fasta')
#print(cline)
clustalw_cline = ClustalwCommandline(clustalw_exe, infile = r'/Users/omidard/Desktop/alignment/seqs.fasta')
assert os.path.isfile(clustalw_exe), "Clustal_W executable is missing"
stdout, stderr = clustalw_cline()
    #print(clustalw_cline)
    #ClustalAlign = AlignIO.read('targetseqs'+str(filenamegenerator[i])+'.aln', "clustal")


#%%print alingment result
#from Bio import AlignIO
ClustalAlign = AlignIO.read(r'/Users/omidard/Desktop/alignment/seqs.aln', "clustal")
print(ClustalAlign) 
tree = Phylo.read(r'/Users/omidard/Desktop/alignment/seqs.dnd', "newick")
Phylo.draw_ascii(tree)
#%%


        
    
#%% des not work 
def find_mutations(recs, ref):
    """Find the mutations in a set of protein records relative to a 
    reference protein sequence"""
    
    mutations = {}
    positions = []
    for rec in recs:
        aln = ClustalAlign = AlignIO.read(r'/Users/omidard/Desktop/alignment/seqs.aln', "clustal")
        #print (aln)
        x = []
        for pos in range(len(aln[0])):
            refaa = aln[0,pos]        
            aa = aln[1,pos]
            if aa != refaa and aa!='-':
                #print (refaa, aln[:,pos], aa)              
                mut = refaa+str(pos+1)+aa
                x.append(mut)
        if len(x)>0:
            mutations[rec.seq] = x
    return mutations

mutations = find_mutations(recs, ref)
mutations.values()


#%%
####pip install biotite
import numpy as np
import biotite
import biotite.sequence as seq
import biotite.sequence.io.fasta as fasta
import biotite.sequence.align as align
import biotite.database.entrez as entrez

# Read each sequence in the file as 'ProteinSequence' object
fasta_file = fasta.FastaFile()
fasta_file.read(r'/Users/omidard/Desktop/alignment/seqs.fasta')
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
import numpy as np
import seaborn as sns
import matplotlib.pylab as plt
  
data_set = identities
ax = sns.heatmap( data_set , linewidth = 0.5 , cmap = 'coolwarm' )
  
plt.title( "2-D Heat Map" )
plt.show()
#%%
# Program to plot 2-D Heat map
# using matplotlib.pyplot.imshow() method
import numpy as np
import matplotlib.pyplot as plt
  
data = identities
plt.imshow( data , cmap = 'autumn' , interpolation = 'nearest' )
  
plt.title( "2-D Heat Map" )
plt.show()
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