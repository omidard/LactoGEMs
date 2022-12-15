#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 17:47:52 2022

@author: omidard
"""

import pandas as pd
from Bio import SeqIO
#%%
def read_value_from_excel(filename, column="C", row=1):

    return pd.read_excel(filename, skiprows=row - 1, usecols=column, nrows=1, header=None, names=["Value"]).iloc[0]["Value"]

#%% run the function to get the gene list for each reaction
num = list(range(1830))
correct_gpr =[]
correct_id =[]
for i in num:
    row = i+3
    react = r'/Users/omidard/Desktop/mar_lbr.xlsx' #excell file of reactome
    reaction_genes = read_value_from_excel (react, 'C',row)
    if ' and ' not in reaction_genes:
        huge = []  
        i1=reaction_genes.replace("(","")
        i2=i1.replace(")","")
        i3=i2.replace('or','\n')
        i4=i3.replace(' ','')
        gene_description = []
        input_genome2 = r'/Users/omidard/Desktop/OTHER/faa/reactome.faa' #merged faa file of 49 lactobacilli genome used for reactome reconstruction
        genome2 = SeqIO.parse(open(input_genome2),'fasta')
        for seq in genome2:
            if seq.id in i4:
                gene_description.append(seq.description)
                short_list = [s.split(None, 1)[1] for s in gene_description]
                short_list2 = list(dict.fromkeys(short_list))
                tgene_collector =[]
                for i in short_list2:
                    gene_name = i #enter the gene name you are looking for
                    collector2 = []
                    input_genome2 = r'/Users/omidard/Desktop/OTHER/faa/reactome.faa'
                    genome2 = SeqIO.parse(open(input_genome2),'fasta')
                    for seq in genome2:
                        if gene_name in seq.description:
                            collector2.append(seq.id)
                    tgene_collector.append(collector2[0])
                huge.append(tgene_collector)
        huge_sh = []
        for i in range(len(huge)):
            for x in huge[i]:
                huge_sh.append(x)
        huge2 = list(dict.fromkeys(huge_sh))
        b = str(huge2)
        b1 = b.replace("['","")
        b2 = b1.replace("']","")
        b3 = b2.replace ("', '"," or ")
        b4 = '(' + b3 + ')'
        correct_gpr.append(b4)
        reaction_ids = read_value_from_excel (react, 'A',row)
        correct_id.append(reaction_ids)
corrected_or = pd.DataFrame({'gpr':correct_gpr,'id':correct_id})
        
#%%
index_collect = []
for i in correct_id:
    index = correct_id.index(i)
    index_collect.append(index)     
        
    
#%%
reaction_genes = pd.read_excel (r'/Users/omidard/Desktop/mar_lbr.xlsx')
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    