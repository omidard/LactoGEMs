#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 27 12:33:53 2022

@author: omidard
"""

import pandas as pd
from Bio import SeqIO
import timeit
#%% run the function to get the gene list for each reaction
start = timeit.default_timer()
#1831
num = list(range(1391))
reaction_genes = pd.read_excel (r'/Users/omidard/Desktop/gpr_reduction/marlbr.xlsx')
correct_gpr =[]
correct_id =[]
for i in num:
    number = str(i)
    if ' and ' not in reaction_genes.GPR[i]:
        reaction_ids = reaction_genes.Abbreviation[i]
        correct_id.append(reaction_ids)
        huge = []  
        i1=reaction_genes.GPR[i].replace("(","")
        i2=i1.replace(")","")
        i3=i2.replace('or','\n')
        i4=i3.replace(' ','')
        gene_description = []
        genome2 = SeqIO.parse(r'/Users/omidard/Desktop/OTHER/faa/reactome.faa','fasta')
        for seq in genome2:
            if seq.id in i4:
                gene_description.append(seq.description)
                short_list = [s.split(None, 1)[1] for s in gene_description]
                short_list2 = list(dict.fromkeys(short_list))
                #print('\n'+'------------------------------','\n')
                #print(short_list2)
                tgene_collector =[]
                for i in short_list2:
                    gene_name = i #enter the gene name you are looking for
                    collector2 = []
                    gleng = []
                    genome2 = SeqIO.parse(r'/Users/omidard/Desktop/OTHER/faa/reactome.faa','fasta')
                    for seq in genome2:
                        if gene_name in seq.description:
                            gleng.append(len(seq))
                            collector2.append(seq.id)
                    len_sorter=pd.DataFrame({'gpr':collector2,'length':gleng})
                    shosh = len_sorter.sort_values(by=['length'],ascending=[False], inplace=False)
                    shoshi = []
                    for i in shosh.gpr:
                        shoshi.append(i)
                    tgene_collector.append(shoshi[0])
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
corrected_or = pd.DataFrame({'gpr':correct_gpr,'id':correct_id})

stop = timeit.default_timer()
print('Time: ', stop - start)  