#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 27 21:53:48 2022

@author: omidard
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 01:17:42 2022

@author: omidard
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 17:47:52 2022

@author: omidard
"""

import pandas as pd
from Bio import SeqIO
import timeit
#%% run the function to get the gene list for each reaction
start = timeit.default_timer()
#1831
num = list(range(1377,1378))
reaction_genes = pd.read_excel (r'/Users/omidard/Desktop/mar_lbr.xlsx')
correct_gpr =[]
correct_id =[]
for i in num:
    if ' and ' in reaction_genes.GPR[i]:
        reaction_ids = reaction_genes.Abbreviation[i]
        correct_id.append(reaction_ids)
        huge = []
        if ')) or' in reaction_genes.GPR[i]:
            i1=reaction_genes.GPR[i].replace(")) or ","\n")
            i2 = i1.splitlines()
            i3 = i2[0].replace(' and ','\n')
            i4 = i3.splitlines() #and containing part
            i5 = i2[1] # or containing part
            
            i4_catch =[]
            for i in range(len(i4)):
                target_querry = i4[i]
                genome2 = SeqIO.parse(r'/Users/omidard/Desktop/multispecies/Supplementary/reactome2.fa','fasta')
                i4_catch2 = []
                for seq in genome2:
                    if seq.id in target_querry:
                        gname = seq.description.split(None, 1)[1]
                        i4_catch2.append(gname)
                i4_catch3 = list(dict.fromkeys(i4_catch2))
                i4_catch.append(i4_catch3)
            for i in range(len(i4_catch)):
                for x in i4_catch[i]:
                    genome2 = SeqIO.parse(r'/Users/omidard/Desktop/multispecies/Supplementary/reactome2.fa','fasta')
                    seqdes = []
                    for seq in genome2:
                        if seq.description in x:
                            seqid = seq.description
                
                
            
#%%                
        i1=reaction_genes.GPR[i].replace(") and (","\n")
        i2=i1.replace("((","")
        i3=i2.replace('))','')
        i4 = i3.splitlines()
        i4_list = []
        for i in range(len(i4)):
            find_list = []
            find = i4[i]
            genome2 = SeqIO.parse(r'/Users/omidard/Desktop/multispecies/Supplementary/reactome2.fa','fasta')
            for seq in genome2:
                if seq.id in find:
                    find_list.append(seq.description)
            i4_list.append(find_list)
        hunter = []    
        for i in range(len(i4_list)):
            x_list = []
            for x in i4_list[i]:
                x2 = x.split(None, 1)[0]
                x_list.append(x2)
            hunter.append(x_list)
            and_list = []
        for i in range(len(hunter)):
            and1 = str(hunter[i])
            and2=and1.replace("', '"," or ")
            and3=and2.replace("['","(")
            and4=and3.replace("']",")")
            and_list.append(and4)
        correct_and = ' and '.join(and_list)
        correct_gpr.append(correct_and)
corrected_and = pd.DataFrame({'gpr':correct_gpr,'id':correct_id})
#%%
          
#%%
        for i in range(len(i4)):
            find = i4[i]
            gene_description = []
            genome2 = SeqIO.parse(r'/Users/omidard/Desktop/multispecies/Supplementary/reactome2.fa','fasta')
            for seq in genome2:
                if seq.id in find:
                    gene_description.append(seq.description)
                    short_list = [s.split(None, 1)[1] for s in gene_description]
                    short_list2 = list(dict.fromkeys(short_list))
                    tgene_collector =[]
                    for i in short_list2:
                        gene_name = i #enter the gene name you are looking for
                        collector2 = []
                        gleng = []
                        genome2 = SeqIO.parse(r'/Users/omidard/Desktop/multispecies/Supplementary/reactome2.fa','fasta')
                        for seq in genome2:
                            if gene_name in seq.description:
                                gleng.append(len(seq))
                                collector2.append(seq.id)
                        len_sorter=pd.DataFrame({'gpr':collector2,'length':gleng})
                        len_sorter.sort_values(by=['length'], inplace=True)
                        tgene_collector.append(len_sorter.gpr[0])
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
#%%
corrected_or.to_excel("output.xlsx")
#%%
index_collect = []
for i in correct_id:
    index = correct_id.index(i)
    index_collect.append(index)     
        
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    