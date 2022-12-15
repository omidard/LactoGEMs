#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 13 23:07:27 2022

@author: omidard
"""
#%%
import pandas as pd
from glob import glob
from Bio import Entrez, SeqIO
import cobra
import os
from os.path import join
from os.path import isfile, join
from os import listdir

#%%
def get_file_names (directory1,directory2):
    onlyfiles = [f for f in listdir(directory1) if isfile(join(directory1, f))]
    gl = str(onlyfiles)
    gl1=gl.replace("', '","\n")
    gl2=gl1.replace(' ','')
    gl3=gl2.replace('[','')
    gl4=gl3.replace(']','')
    gl5=gl4.replace("'","")
    gl6=gl5.replace(".gbk","")
    gl7=gl6.replace('.DS_Store','')
    with open(directory2+'/targ.txt','w') as tf:
        tf.writelines(gl7)
        tf.close()
    gl8 = []
    txtf = directory2+'/targ.txt'
    with open(txtf) as f:
        for line in f:
            line = line.strip()
            if line != "":
                gl8.append(line)
    gl9=pd.DataFrame({'models':gl8})
    path = directory2+'/modelInformation.xlsx'
    modelInformation = gl9.to_excel(path)
#%%
directory1 =r'/Users/omidard/Desktop/MULTI/Models'
directory2 =r'/Users/omidard/Desktop/multispecies/Supplementary'
get_file_names (directory1,directory2)
#%%
modelOfInterest=pd.read_excel(directory2+'/modelInformation.xlsx')
modelOfInterest.pop('Unnamed: 0')
print(modelOfInterest)
#%%
# gather the general information on the draft models
modelid = []
genumb = []
reanumb = []
for strain in modelOfInterest.models:
    model=cobra.io.load_json_model(directory1+'/'+strain)
    print (model.id,'Number of Model Genes:',len(model.genes),'Number of Model Reactions:',len(model.reactions))
    modelid.append(model.id)
    genumb.append(len(model.genes))
    reanumb.append(len(model.reactions))
modelsinf = pd.DataFrame({'model_id':modelid,'reactions':reanumb,'genes':genumb})
#%%
directory1 =r'/Users/omidard/Desktop/multispecies/Supplementary/genomes'
files = glob('%s/*.gbk'%directory1)
strain_info = []
for file in files:
        handle = open(file)
        record = list(SeqIO.parse(handle, "genbank"))
        for i in record:
            for f in i.features:
                if f.type=='source':
                    info = {}
                    info['file'] = file.replace('/Users/omidard/Desktop/multispecies/Supplementary/genomes/','')
                    info['id'] = file.replace('/Users/omidard/Desktop/multispecies/Supplementary/genomes/','')
                    for q in f.qualifiers.keys():
                        info[q] = '|'.join(f.qualifiers[q])
                        strain_info.append(info)
                        sinf = pd.DataFrame(strain_info)
sinf2 = sinf.drop_duplicates(subset='id', keep='first', inplace=False, ignore_index=False)
#%%
organism = list(sinf2.organism)
strain = list(sinf2.strain)
modelsinf = pd.DataFrame({'organism':organism,'strain':strain,'model_id':modelid,'reactions':reanumb,'genes':genumb})
modelsinf.pop('model_id')
modelsinf.sort_values(by=['reactions'])
print(modelsinf.sort_values(by=['reactions']))
#%%
import matplotlib.pyplot as plt
modinf = modelsinf.sort_values(by=['reactions'])
print(modinf)
plt.figure(figsize=(100, 50))
plt.ylabel('Number of reactions',fontsize=100, color='red')
plt.xlabel('Models',fontsize=100, color='red')
plt.rc('xtick', labelsize=70) 
plt.rc('ytick', labelsize=70)
plt.xticks(range(len(modinf.organism)), modinf.organism, rotation='vertical')
plt.margins(0.05)
plt.subplots_adjust(bottom=0.15)
plt.xticks(fontsize= 25)

plt.grid(linestyle="--", linewidth=0.95, color='.25', zorder=-10)
plt.bar(modinf.organism, modinf.reactions)
print (modinf.organism[17])
print(len(modinf.organism))
print(len(modinf.reactions))
#%%target genome information

directory1 =r'/Users/omidard/Desktop/MULTI/multispecies3/Supplementary/genomes'
files = glob('%s/*.gbk'%directory1)
strain_info = []
for file in files:
        handle = open(file)
        record = list(SeqIO.parse(handle, "genbank"))
        for i in record:
            for f in i.features:
                if f.type=='source':
                    info = {}
                    info['file'] = file.replace('/Users/omidard/Desktop/multispecies/Supplementary/genomes/','')
                    info['id'] = file.replace('/Users/omidard/Desktop/multispecies/Supplementary/genomes/','')
                    for q in f.qualifiers.keys():
                        info[q] = '|'.join(f.qualifiers[q])
                        strain_info.append(info)
                        sinf = pd.DataFrame(strain_info)
sinf2 = sinf.drop_duplicates(subset='id', keep='first', inplace=False, ignore_index=False)
print(sinf2)
#%%
directory1 =r'/Users/omidard/Desktop/multispecies/Supplementary/genomes'
files = glob('%s/*.gbk'%directory1)
strain_info = []
for file in files:
        handle = open(file)
        record = list(SeqIO.parse(handle, "genbank"))

#%%
result_file = r'/Users/omidard/Desktop/multispecies/Supplementary/prots/GCF1_009556455.1.fa'
fasta_file = r'/Users/omidard/Desktop/multispecies/Supplementary/prots/GCF_009556455.1.fa'
with open(result_file, "w") as g:
    fasta_sequences = SeqIO.parse(fasta_file,'fasta')
    for seq in fasta_sequences:
            SeqIO.write([seq], g, "fasta")
            
            
#%%
model = cobra.io.load_matlab_model(os.path.join("LBR3.mat"))
#%%
## Load the previously generated homology matrix for E. coli strains of interest
hom_matrix=pd.read_csv('ortho_matrix.csv')
hom_matrix=hom_matrix.set_index('Unnamed: 0')
#%%
import cobra
import pandas as pd
from cobra.io import load_json_model
from glob import glob
from cobra.manipulation.delete import delete_model_genes, remove_genes
import cobra
import os
from os.path import join
#%%
#create strain-specific draft models and save them
for strain in hom_matrix.columns:
    
    #Get the list of Gene IDs from the homology matrix dataframe for the current strain without a homolog
    currentStrain=hom_matrix[strain]
    nonHomologous=currentStrain[currentStrain==0.0]
    nonHomologous=nonHomologous.index.tolist()
    nonHomologous.remove('spontaneous')
    nonHomologous.remove('EXCHANGE')
    nonHomologous.remove('BIOMASS')
    nonHomologous.remove('Diffusion')
    nonHomologous.remove('GAP')
    nonHomologous.remove('DEMAND')
    #above genes are artificial genes used in reactome for spontaneous/exchange/gap/biomass/demand/diffusion reactions and as such has no homologs,
    #However, it is retained for these spontaneous reactions to function
    #nonHomologous.remove('s0001')
    #Define a list of Gene objects from the base reconstruction to be deleted from the current strain
    toDelete=[]
    for gene in nonHomologous:
        toDelete.append(model.genes.get_by_id(gene))

    #Establish a model copy and use the COBRApy function to remove the appropriate content and save this model
    modelCopy=model.copy()
    remove_genes(modelCopy, toDelete, remove_reactions=True)
    modelCopy.id=str(strain)
    cobra.io.json.save_json_model(modelCopy, str('Models/'+strain+'.json'), pretty=False)        
#%%
models=glob('%s/*.json'%'Models')
geneIDs_matrix=pd.read_csv('geneIDs_matrix.csv')
geneIDs_matrix=geneIDs_matrix.set_index('Unnamed: 0')
geneIDs_matrix
#%%
from cobra.manipulation.modify import rename_genes

for mod in models:
    model=cobra.io.load_json_model(mod)
    for column in geneIDs_matrix.columns:
        if column in mod.replace('.json',''):
            currentStrain=column
    
    IDMapping=geneIDs_matrix[currentStrain].to_dict()
    IDMappingParsed = {k:v for k,v in IDMapping.items() if v != 'None'}
    
    rename_genes(model,IDMappingParsed)
    cobra.io.json.save_json_model(model,mod, pretty=False)




