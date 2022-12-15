#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 00:11:53 2022

@author: omidard
"""

import pandas as pd
from Bio import SeqIO

#%%function for retriving locus tags for specific reactions in reactome


def read_value_from_excel(filename, column="C", row=1):

    return pd.read_excel(filename, skiprows=row - 1, usecols=column, nrows=1, header=None, names=["Value"]).iloc[0]["Value"]


#%% run the function to get the gene list for each reaction
row = 15 #change the numerical value based on reaction's row in excell


react = r'/Users/omidard/Desktop/mar_lbr.xlsx' #excell file of reactome
reaction_genes = read_value_from_excel (react, 'C',row)
and_seprator = reaction_genes.replace('and','\n')
and_seprator2 = and_seprator.replace(' (','(')
and_seprator3 = and_seprator2.replace('((','(')
and_seprator4= and_seprator3.replace('))',')')
#and_seprator5=and_seprator4.replace(' ','')
with open ('/Users/omidard/Desktop/gpr_reduction/gprs.txt','w') as w:
    w.writelines(and_seprator4)
with open ('/Users/omidard/Desktop/gpr_reduction/gprs.txt','r') as r:
    gprs = set()
    for line in r:
        line = line.strip()
        if line != "":
            gprs.add(line)
huge = []           
gprs2 = list(gprs)
for z in gprs2:
        i1=z.replace("(","")
        i2=i1.replace(")","")
        i3=i2.replace('or','\n')
        i4=i3.replace(' ','')
        gene_description = []
        input_genome2 = r'/Users/omidard/Desktop/OTHER/faa/reactome.faa' #merged faa file of 49 lactobacilli genome used for reactome reconstruction
        genome2 = SeqIO.parse(open(input_genome2),'fasta')
        for seq in genome2:
            if seq.id in i4:
                #print(seq)
                #print(seq.description)
                gene_description.append(seq.description)
        short_list = [s.split(None, 1)[1] for s in gene_description]
        short_list2 = list(dict.fromkeys(short_list))
        gprs3 =[]
        tgene_collector =[]
        for i in short_list2:
            gene_name = i #enter the gene name you are looking for
            collector2 = []
            input_genome2 = r'/Users/omidard/Desktop/OTHER/faa/reactome.faa'
            genome2 = SeqIO.parse(open(input_genome2),'fasta')
            for seq in genome2:
                 #if gene_name in seq.description:# and 'subunit' not in seq.description:
                 if gene_name in seq.description:# and 'Phenyl' not in seq.description:# and 'D-' not in seq.description and 'Methionine' not in seq.description:
                    #print(seq.id)
                    #print(seq.description)
                    collector2.append(seq.id)
                    #print(collector2)
            tgene_collector.append(collector2)
        huge.append(tgene_collector)
bbb=[]
for i in range(len(huge)):
    fof = huge[i]
    #print(fof)
    aaa = []
    for x in range(len(fof)):
        fof2 = fof[x][0]
        #print(fof2)
        gepre = str(fof2)
        aaa.append(gepre)
    bbb.append(aaa)

b3 = []
for i in range(len(bbb)):
    b2 = str(bbb[i])
    b3.append(b2)

b8 = []
for i in range(len(b3)):
    b4 = b3[i].replace("['","")
    b5 = b4.replace("']","")
    b6 = b5.replace ("', '"," or ")
    b7 = '(' + b6 + ')'
    b8.append(b7)
b9 = ' and '.join(b8)
#%%       
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
#%%                       
        gene_and_locus = pd.DataFrame({'gene_names':to_print,'locus_tags':collector2})
        gprs3.append(gene_and_locus)
        tgene = gprs3[0].locus_tags[0]
        tgene_collector.append(tgene)

#%%
print('\n','\n','-----------------------------------------------------------------','\n','full name of genes involved in the targeted reaction','\n','-----------------------------------------------------------------','\n')
#% find gene ids for manual boolean rule curation
#reaction_genes = 'AADJCKDE_00878','ADJCKDE_01855','AADJCKDE_01854'
gene_description = []
involved_genes =[]
input_genome2 = r'/Users/omidard/Desktop/OTHER/faa/reactome.faa' #merged faa file of 49 lactobacilli genome used for reactome reconstruction
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


#print('\n','#--------------------------------------------------------------#')
#print(' #--------------------reaction information----------------------#')
#print(' #--------------------------------------------------------------#')
#print('\n','>>>>>>>>>','Modelseed_id =',read_value_from_excel (react, 'O', row))
#print('\n','>>>>>>>>>','Reaction =',read_value_from_excel (react, 'B', row))
#print('\n','>>>>>>>>>',"Reaction's name = ",read_value_from_excel (react, 'D', row))
#print('\n','-----------------------------------------------------------------','\n','>>>>> different type of genes has been found in this reaction are <<<<<<<','\n','-----------------------------------------------------------------','\n',statement4,'\n','-----------------------------------------------------------------','\n')
#print ('\n','!!!!!! next step instruction !!!!!!','\n','------>>> in cell below change the value for <T> variable from [[[ 0 to',len(involved_genes)-1,']]]<<<--------','\n','\n','\n','-----------------------------------------------------------------','\n','\n')
#%%
#gprs = []
for i in involved_genes:
    gene_name = i #enter the gene name you are looking for
    involved_genes3 = []
    gene_description2 =[]
    collector2 = []
    input_genome2 = r'/Users/omidard/Desktop/OTHER/faa/reactome.faa'
    genome2 = SeqIO.parse(open(input_genome2),'fasta')
    for seq in genome2:
         #if gene_name in seq.description:# and 'subunit' not in seq.description:
         if gene_name.lower() in seq.description.lower():# and 'Phenyl' not in seq.description:# and 'D-' not in seq.description and 'Methionine' not in seq.description:
            #print(seq.id)
            print(seq.description)
            collector2.append(seq.id)
            gene_description2.append(seq.description)
    TO_PRINT = [s.split(None, 1)[1] for s in gene_description2]
    gene_and_locus = pd.DataFrame({'gene_names':TO_PRINT,'locus_tags':collector2})
    gprs.append(gene_and_locus)
    TO_PRINT2 = pd.DataFrame(TO_PRINT)
    TO_PRINT2.columns = ['enzymes']
    TO_PRINT3 = list(TO_PRINT2.groupby('enzymes'))
    print(TO_PRINT3)
    
#%%
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
    print('\n',search_result_v3,'\n','=========================================','\n','\n','\n')

#%%
def get_gene_lens(query, in_folder='prots'):

    file = '%s/%s.faa'%(in_folder, query)
    handle = open(file)
    records = SeqIO.parse(handle, "fasta")
    out = []
    
    for record in records:
        out.append({'gene':record.name, 'gene_length':len(record.seq)})
    
    out = pd.DataFrame(out)
    return out
#%%
genes_length = get_gene_lens('reactome', in_folder=r'/Users/omidard/Desktop/OTHER/faa')
#%% collect similar genes from lactobacilli metagenome
inf_colc = []
for i in range(len(involved_genes)):
    target = involved_genes[i]
    print(target)
    collector = []
    input_genome2 = r'/Users/omidard/Desktop/other/faa/reactome.faa'
    genome2 = SeqIO.parse(open(input_genome2),'fasta')
    for seq in genome2:
        if target in seq.description: #and 'YjjG' not in seq.description: 
            #print(seq.id)
            #print(seq.description)
            collector.append(seq.id)
            g_l=[]
            for i in genes_length:
                g_l.append(genes_length[i])
            ge = []
            gl = []       
            for g in collector:
                if str(g) in list(g_l[0]):
                    index = list(g_l[0]).index(g)
                    ge.append(g)
                    gl.append(g_l[1][index])
                    gene_inf = pd.DataFrame({'gene':ge,'gene_length':gl})
                    inf_colc.append(gene_inf)
    can_colc = []
    for i in range(len(inf_colc)):
        liss = list(inf_colc[i].groupby('gene_length'))
        candid_gene = liss[len(liss)-1][1].gene
        can_colc.append(candid_gene)
        
    #print('\n','---------------->>>>>>','\n','number of found genes in genome = ',len(collector),'\n','number of genes in reactome =',reaction_genes.count('_'), '<<<<----------------')
    gene_set = str(collector)
    target_locustag_v1 = gene_set.replace("', '"," or ")
    target_locustag_v2 = target_locustag_v1.replace("['","(")
    target_locustag_v3 = target_locustag_v2.replace("']",")")
    #print('\n','-----------------------------------------------------------------','\n','-----------------------------------------------------------------')
    #print('\n',target_locustag_v3,'\n','---------------------------------------------')

#%%
#__________________searching for a gene name in reactome FASTA file---------

gene_name = "ketolase" #enter the gene name you are looking for
involved_genes3 = []
gene_description2 =[]
collector2 = []
input_genome2 = r'/Users/omidard/Desktop/faa/reactome.faa'
genome2 = SeqIO.parse(open(input_genome2),'fasta')
for seq in genome2:
     #if gene_name in seq.description:# and 'subunit' not in seq.description:
     if gene_name.lower() in seq.description.lower():# and 'Phenyl' not in seq.description:# and 'D-' not in seq.description and 'Methionine' not in seq.description:
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

#%%
def get_gene_lens(query, in_folder='prots'):

    file = '%s/%s.faa'%(in_folder, query)
    handle = open(file)
    records = SeqIO.parse(handle, "fasta")
    out = []
    
    for record in records:
        out.append({'gene':record.name, 'gene_length':len(record.seq)})
    
    out = pd.DataFrame(out)
    return out

