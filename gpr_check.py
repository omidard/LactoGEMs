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
row = 767 #change the numerical value based on reaction's row in excell


react = r'/Users/omidard/Desktop/mar_lbr_apr.xlsx'
#react = r'/Users/omidard/Desktop/gpr_reduction/marlbr.xlsx'
reaction_genes = read_value_from_excel (react, 'C',row)
#print('\n','\n','-----------------------------------------------------------------','\n','full name of genes involved in the targeted reaction','\n','-----------------------------------------------------------------','\n')
#% find gene ids for manual boolean rule curation
#reaction_genes = 'AADJCKDE_00878','ADJCKDE_01855','AADJCKDE_01854'
gene_description = []
involved_genes =[]
input_genome2 = r'/Users/omidard/Desktop/OTHER/faa/reactome.fa' #merged faa file of 49 lactobacilli genome used for reactome reconstruction
genome2 = SeqIO.parse(open(input_genome2),'fasta')
for seq in genome2:
    if seq.id in reaction_genes:
        print(len(seq))
        print(seq)
        #print(seq.description)
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


print('\n','#--------------------------------------------------------------#')
print(' #--------------------reaction information----------------------#')
print(' #--------------------------------------------------------------#')
print('\n','>>>>>>>>>','Modelseed_id =',read_value_from_excel (react, 'O', row))
print('\n','>>>>>>>>>','Reaction =',read_value_from_excel (react, 'B', row))
print('\n','>>>>>>>>>',"Reaction's name = ",read_value_from_excel (react, 'D', row))
print('\n','-----------------------------------------------------------------','\n','>>>>> different type of genes has been found in this reaction are <<<<<<<','\n','-----------------------------------------------------------------','\n',statement4,'\n','-----------------------------------------------------------------','\n')
print ('\n','!!!!!! next step instruction !!!!!!','\n','------>>> in cell below change the value for <T> variable from [[[ 0 to',len(involved_genes)-1,']]]<<<--------','\n','\n','\n','-----------------------------------------------------------------','\n','\n')

#%% collect similar genes from lactobacilli metagenome
T = 0
target = involved_genes[T] #change the value based on printed instruction to collect involved locus tags in this reaction for all possible alleles of gene of interest
collector = []
input_genome2 = r'/Users/omidard/Desktop/OTHER/faa/reactome.fa'
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

gene_name = "phosphoadenosine" #enter the gene name you are looking for
involved_genes3 = []
gene_description2 =[]
collector2 = []
input_genome2 = r'/Users/omidard/Desktop/OTHER/faa/reactome.fa'
genome2 = SeqIO.parse(open(input_genome2),'fasta')
for seq in genome2:
     #if gene_name in seq.description:# and 'subunit' not in seq.description:
     if gene_name.lower() in seq.description.lower():# and 'Phenyl' not in seq.description:# and 'D-' not in seq.description and 'Methionine' not in seq.description:
        #print(seq.id)
        print(seq.description)
        print(len(seq))
        collector2.append(seq.id)
        gene_description2.append(seq.description)
TO_PRINT = [s.split(None, 1)[1] for s in gene_description2]
TO_PRINT2 = pd.DataFrame(TO_PRINT)
TO_PRINT2.columns = ['enzymes']
TO_PRINT3 = list(TO_PRINT2.groupby('enzymes'))

for i in range(len(TO_PRINT3)):
    desc2 = TO_PRINT3[i][0]
    involved_genes3.append(desc2)
statementv1 = str(involved_genes3)
statementv2 = statementv1.replace("['",'\n')
statementv3 = statementv2.replace("', '","\n")
statementv4 = statementv3.replace("']","\n")


print('\n','-------------------------------------------------','\n','>>>>gene types containing search keyword<<<<','\n',statementv4)
print('\n','----------- number of found genes = ',len(collector2),'----------------')
gene_set2 = str(collector2)
search_result_v1 = gene_set2.replace("', '"," or ")
search_result_v2 = search_result_v1.replace("['","(")
search_result_v3 = search_result_v2.replace("']",")")
print('\n','!!!!!!!instruction!!!!!!!','\n','copy the below locustags list and paste it in the related cell in reactome excell shit ')
print('\n',search_result_v3)











#%%
#other things

#%%
#matching locus tags for leru and reactome
leru_genome = r'/Users/omidard/Desktop/multispecies/Supplementary/prots/lreu.fa'

lreugenome = SeqIO.parse(open(leru_genome),'fasta')
for seq in lreugenome:
    print(seq.description)
    
#%% =============unassigned finder================   
#%% unassigned metabolic genes finder: preparing a list of genes in reactome

collector =[]
reactome_genome = r'/Users/omidard/Desktop/faa/reactome.faa'
reactome = pd.read_excel(r'/Users/omidard/Desktop/7-jan-LbReactome2_2.xlsx')
tags = reactome.GPR[:90000000]
for i in tags:
    tags2 = str(i)[:1000000]
    collector.append(tags2)
tags3=str(collector)
lis1=tags3.replace("', '","\n")
lis2=lis1.replace("'","")
lis3=lis2.replace('(','')
lis4=lis3.replace(')','')
lis5=lis4.replace('and','\n')
lis6=lis5.replace('or','\n')
lis7=lis6.replace('\ufeff','')
lis8=lis7.replace(' ','')
lis9=lis8.replace('[','')
lis10=lis9.replace(']','')


#%% find all unassigned genes from reactome fasta
unassigned=[]
reactome_genome2 = SeqIO.parse(open(reactome_genome),'fasta')  
for seq in reactome_genome2:
    #print(seq.id)
    if seq.id not in lis10:
        unassigned.append(seq.description)

flis = [s.split(None, 1)[1] for s in unassigned]
flis2 = [s.split(None, 1)[0] for s in unassigned]
flis3= pd.DataFrame(flis)
flis3.columns = ['enzymes']
flis4=list(flis3.groupby('enzymes'))
final_list=flis4.sort()
#%% removing non-metabolic genes
targets=[]
for i in range(len(flis4)):
    seed=flis4[i][0]
    if 'transposase' not in seed.lower() and 'dna ligase' not in seed.lower() and 'chaperonin' not in seed.lower() and 'heat shock' not in seed.lower() and 'ribosomal protein' not in seed.lower() and 'helicase' not in seed.lower() and 'CRISPR' not in seed and 'repressor' not in seed and 'regulator' not in seed and 'RNA' not in seed and 'DNA' not in seed and 'Flagel' not in seed and 'stress' not in seed and 'Regulatory' not in seed and 'Transcription' not in seed and 'Signal' not in seed and 'ribosomal' not in seed and 'estriction' not in seed and 'multidrug' not in seed.lower() and 'Accessory' not in seed and 'Inner membrane protein' not in seed and 'Cell division' not in seed and 'catabolite' not in seed.lower() and 'Cold shock' not in seed and 'Elongation' not in seed and 'Bac' not in seed and 'activator' not in seed and 'Protein' not in seed and 'Trans' not in seed and 'Endo' not in seed and 'Ribonuclease' not in seed and 'resistance' not in seed and 'Nucleoid' not in seed and 'bosom' not in seed and 'rotease' not in seed and 'ensor' not in seed and 'receptor' not in seed and 'eplication' not in seed:
        targets.append(seed)
#%%
number=958
print('\n','\n','\n','\n',targets[number],'\n','\n','\n','\n')

#%%search for previuosly unassigned genes locus tags
gene_name = targets[number] #enter the gene name you are looking for
involved_genes3 = []
gene_description2 =[]
collector2 = []
input_genome2 = r'/Users/omidard/Desktop/faa/reactome.faa'
genome2 = SeqIO.parse(open(input_genome2),'fasta')
for seq in genome2:
     #if gene_name in seq.description:# and 'subunit' not in seq.description:
     if gene_name.lower() in seq.description.lower() and 'Bifunctional' not in seq.description:# and 'D-' not in seq.description and 'Methionine' not in seq.description:
        #print(seq.id)
        print(seq.description)
        collector2.append(seq.id)
        gene_description2.append(seq.description)
TO_PRINT = [s.split(None, 1)[1] for s in gene_description2]
TO_PRINT2 = pd.DataFrame(TO_PRINT)
TO_PRINT2.columns = ['enzymes']
TO_PRINT3 = list(TO_PRINT2.groupby('enzymes'))

for i in range(len(TO_PRINT3)):
    desc2 = TO_PRINT3[i][0]
    involved_genes3.append(desc2)
statementv1 = str(involved_genes3)
statementv2 = statementv1.replace("['",'\n')
statementv3 = statementv2.replace("', '","\n")
statementv4 = statementv3.replace("']","\n")


print('\n','-------------------------------------------------','\n','>>>>gene types containing search keyword<<<<','\n',statementv4)
print('\n','----------- number of found genes = ',len(collector2),'----------------')
gene_set2 = str(collector2)
search_result_v1 = gene_set2.replace("', '"," or ")
search_result_v2 = search_result_v1.replace("['","(")
search_result_v3 = search_result_v2.replace("']",")")
print('\n','!!!!!!!instruction!!!!!!!','\n','copy the below locustags list and paste it in the related cell in reactome excell shit ')
print('\n',search_result_v3)

#%% extract unassigned genes to a new fasta file
input_genome2 = r'/Users/omidard/Desktop/faa/reactome.faa' #merged faa file of 49 lactobacilli genome used for reactome reconstruction
genome2 = SeqIO.parse(open(input_genome2),'fasta')
f = open('unassigned.fasta','w')
for seq in genome2:
    if seq.id in flis2:
        SeqIO.write([seq], f, "fasta")

#%% extract locus tags based on a list of gene description produced by blast search result
lreu_description=[]
lreu_genome = r'/Users/omidard/Desktop/LBpubgems/lreu.faa' #merged faa file of 49 lactobacilli genome used for reactome reconstruction
lreu_genome2 = SeqIO.parse(open(lreu_genome),'fasta')
glist = r'/Users/omidard/Desktop/LBpubgems/teglist2.xlsx'
gok=[]
glist2= pd.read_excel(glist)
for i in glist2.lreu_genes:
    glist3= str(i)[:1000000]
    gok.append(glist3)
gok1=str(gok)
gok2=gok1.replace("', '","\n")
gok3=gok2.replace("['","")
gok4=gok3.replace("']","")
for seq in lreu_genome2:
    query=seq.description.split(None, 1)[0]
    if query in gok:
        lreu_description.append(seq.description)
#%%extract exact locus tags    
locus_tags=[]  
for i in lreu_description:
    fis1= i.split(None, 1)[1]
    fis2=fis1.split(None, 1)[0]
    locus_tags.append(fis2)
loc=str(locus_tags)
loc2=loc.replace("]', '[","\n")
loc3=loc2.replace("['[","")
loc4=loc3.replace("]']","")
loc5=loc4.replace('locus_tag=','') #locus tags list
#%%FIND REACTIONS which has specific locus tags
import numpy as np
newloc = r'/Users/omidard/Desktop/LBpubgems/newloc.xlsx'
model = r'/Users/omidard/Desktop/LBpubgems/Additional_file_6d.xlsx'
fook=[]
newloc2= pd.read_excel(newloc)
for i in newloc2.LTAG:
    fok= str(i)[:1000000]
    fook.append(fok)
fook1=str(fook)
fook2=fook1.replace("', '","\n")
fook3=fook2.replace("['","")
fook4=fook3.replace("']","") #final tag list
model2= pd.read_excel(model)
model3=model2.GPR
newloc3=newloc2.LTAG
for i in range(len(newloc2)): 
    if newloc3[i] in str(model3):
        indx= np.where(model3 == newloc3[i])
        print(model2.ID[indx])
        
        


