# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
#import required mudules
#pip install funcy
import pandas as pd
import os
import numpy as np
import pickle
#%%
f = open('planti_gene_map_rev.txt', 'rb')
genmap = pickle.load(f)
#dividing the genmap dictionary into keys and values list
gmkeys = list(genmap.keys())
gmval = list(genmap.values())
#load excel file
gpr = pd.read_excel ('ugpr.xls')
print (gpr)
#%%turn ugpr dataframe to a list
gprlist = gpr.GPR.tolist()
#%%
#turn a list into a single string
listToStr = ' '.join([str(elem) for elem in gprlist])
#TURN STrINGS WITH SPACE INTO A LIST
gpr3 = listToStr.split()
#%% reduce dict
from funcy import project
redict = project(genmap, gpr3)
query = list(redict.values())
#%%grouping dict based on similar values
from collections import defaultdict
genmapgrouped = defaultdict(list)
for key, value in sorted(genmap.items()):
    genmapgrouped[value].append(key)
#%%reduce dict
groupedgenes = project(genmapgrouped, query)
for keys in groupedgenes:
    sep = str(groupedgenes[keys]).replace(' ', '\r\n')
    with open (str(keys)+'.txt',"w") as big:
        sep2 = sep.replace(",", "")
        sep3 = sep2.replace("[","")
        sep4 = sep3.replace("]","")
        sep5 = sep4.replace("'","")
        big.write(sep5)
        big.close()

##############################END END END END#################################



#%% alignment input preparation
from Bio import SeqIO
fasta_file = 'pan.faa' # Input fasta file
filenamegenerator2 = list(groupedgenes.keys())
joker2 = list(range(len(list(groupedgenes.keys()))))               
for i in joker2:
    wanted_file2 = str(filenamegenerator2[i])+'.txt'
    result_file2 = 'ver6'+str(filenamegenerator2[i])+'.fasta'
    wanted2 = set()
    with open(wanted_file2,"r") as fati:
        for line in fati:
            line = line.strip()
            if line !="":
                wanted2.add(line)
                with open(result_file2, "w") as g:
                    fasta_sequences2 = SeqIO.parse(open(fasta_file),'fasta') 
                    for seq in fasta_sequences2:
                        if seq.id in wanted2:
                            SeqIO.write([seq], g, "fasta")
                            
#%%step6: run MSA
filenamegenerator2 = list(groupedgenes.keys())
joker2 = list(range(len(list(groupedgenes.keys()))))                                 
fasta_file = "pan.faa" # Input fasta file
for i in joker2:
    wanted_file = str(filenamegenerator2[i])+'.txt'
    result_file = 'ver6'+str(filenamegenerator2[i])+'.fasta'
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
                            
                            
                            
#%%
wanted_file = "C8840.txt"
result_file = "ver6C8840.fasta"   
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
                        
                        
#%%
clustalw_exe = "/Users/omidard/opt/anaconda3/envs/clustalw2/bin/clustalw2"

from Bio.Align.Applications import ClustalwCommandline
cline = ClustalwCommandline("clustalw2", infile = "ver6C10651.fasta")
   
clustalw_cline = ClustalwCommandline(clustalw_exe, infile ="ver6C10651.fasta" )
assert os.path.isfile(clustalw_exe), "Clustal_W executable is missing"
stdout, stderr = clustalw_cline()
#%%get values(cdhits) for a list of keys(locustags)
gennames = [genmap[x] for x in intersection_as_list]
#%%get keys(locustags) for values(cdhits)
def getKeysByValue(dictOfElements, valueToFind):
    listOfKeys = list()
    listOfItems = dictOfElements.items()
    for item  in listOfItems:
        if item[1] == valueToFind:
            listOfKeys.append(item[0])
    return  listOfKeys
numbank =  list(range(len (gennames)))
for i in numbank:
    listOfKeys = getKeysByValue(genmap, gennames[i])
    for key  in listOfKeys:
        with open ('cdhitsmapingresults.txt', "a") as m:
            m.writelines(str(listOfKeys))
    with open ('cdhitsmapingresults.txt', "r") as g:
     alllactogenes = list(g)
     
     
#%%
from Bio import SeqIO

with open('pan.faa','r') as fasta_file:
    record_dict = SeqIO.to_dict(open(SeqIO.parse(fasta_file, 'fasta')))

#with open('C31413.txt','r') as text_file:
   # orthologs_txt = text_file.read()

#genes_to_keep = []     
for ortholog in orthologs_txt.splitlines():
    try:
        genes_to_keep.append( record_dict[ortholog] )
    except KeyError:
        pass

with open('ortholog1.fasta','w') as output_file:
    SeqIO.write(genes_to_keep, output_file, 'fasta')
#%%
#Find the indices at which any element of one list occurs in another
st = set(intersection_as_list)
indices_list = [i for i, e in enumerate(cl) if e in st]
#extract cdhits using mapped indices
cdhits = list(genmap.items())


