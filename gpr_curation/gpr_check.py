#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 27 10:02:14 2022

@author: omidard
"""
'''script for semi-manual curation of boolean rules in reactome'''

#%%
import pandas as pd
from Bio import SeqIO

#%%function for retriving locus tags for specific reactions in reactome


def read_value_from_excel(filename, column="C", row=1):

    return pd.read_excel(filename, skiprows=row - 1, usecols=column, nrows=1, header=None, names=["Value"]).iloc[0]["Value"]

print('\n','\n','-----------------------------------------------------------------','\n','!!!!!!!instruction!!!!!!!','\n','change the numerical value in below code based on target reactions row in excell shit','\n','-----------------------------------------------------------------','\n','\n','\n','\n',)


#%% run the function to get the gene list for each reaction
row = 83 #change the numerical value based on reaction's row in excell
react = r'/Users/omidard/Desktop/7-jan-LbReactome2_2.xlsx' #excell file of reactome
reaction_genes = read_value_from_excel (react, 'C',row ) 
print('\n','\n','-----------------------------------------------------------------','\n','full name of genes involved in the targeted reaction','\n','-----------------------------------------------------------------','\n')
#% find gene ids for manual boolean rule curation
#reaction_genes = 'AADJCKDE_00878','ADJCKDE_01855','AADJCKDE_01854'
gene_description = []
involved_genes =[]
input_genome2 = r'/Users/omidard/Desktop/faa/reactome.faa' #merged faa file of 49 lactobacilli genome used for reactome reconstruction
genome2 = SeqIO.parse(open(input_genome2),'fasta')
for seq in genome2:
    if seq.id in reaction_genes:
        #print(seq)
        print(seq.description)
        gene_description.append(seq.description)
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
print('\n','------------------------------------------------------------')
print('------------------reaction information-----------------------')
print('\n','>>>>>>>>>','Modelseed_id =',read_value_from_excel (react, 'O', row))
print('\n','>>>>>>>>>','Reaction =',read_value_from_excel (react, 'B', row))
print('\n','>>>>>>>>>',"Reaction's name = ",read_value_from_excel (react, 'D', row))
print('\n','-----------------------------------------------------------------','\n','>>>>> different type of genes has been found in this reaction are <<<<<<<','\n',statement4,'\n','-----------------------------------------------------------------','\n')
print ('\n','!!!!!! next step instruction !!!!!!','\n','------>>> in cell below change the value in <target> variable from [[[ 0 to',len(involved_genes)-1,']]]<<<--------','\n','\n','\n','-----------------------------------------------------------------','\n','\n')



#%% collect similar genes from lactobacilli metagenome

target = involved_genes[0] #change the value based on printed instruction to collect involved locus tags in this reaction for all possible alleles of gene of interest
collector = []
input_genome2 = r'/Users/omidard/Desktop/faa/reactome.faa'
genome2 = SeqIO.parse(open(input_genome2),'fasta')
for seq in genome2:
    if target in seq.description: #and 'YjjG' not in seq.description: 
        #print(seq.id)
        print(seq.description)
        collector.append(seq.id)
print('\n','---------------->>>>>> number of found genes in genome = ',len(collector),'\n','number of genes in reactome =',reaction_genes.count('_'), '<<<<----------------')
gene_set = str(collector)
target_locustag_v1 = gene_set.replace("', '"," or ")
target_locustag_v2 = target_locustag_v1.replace("['","(")
target_locustag_v3 = target_locustag_v2.replace("']",")")
print('\n','!!!!!!!instruction!!!!!!!','\n','copy the below locustags list and paste it in the related cell in reactome excell shit ')
print('\n',target_locustag_v3)



#%%

#__________________searching for a gene name in reactome FASTA file---------

gene_name =  "malic" #enter the gene name you are looking for
collector2 = []
input_genome2 = r'/Users/omidard/Desktop/faa/reactome.faa'
genome2 = SeqIO.parse(open(input_genome2),'fasta')
for seq in genome2:
    if gene_name in seq.description: #and 'YjjG' not in seq.description:
        #print(seq.id)
        print(seq.description)
        collector2.append(seq.id)
print('\n','----------- number of found genes = ',len(collector2),'----------------')
gene_set2 = str(collector2)
search_result_v1 = gene_set2.replace("', '"," or ")
search_result_v2 = search_result_v1.replace("['","(")
search_result_v3 = search_result_v2.replace("']",")")
print('\n','!!!!!!!instruction!!!!!!!','\n','copy the below locustags list and paste it in the related cell in reactome excell shit ')
print('\n',search_result_v3)










