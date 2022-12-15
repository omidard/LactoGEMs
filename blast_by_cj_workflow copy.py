#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 11 00:08:40 2022

@author: omidard

"""

#%%
#import packages needed
from os import listdir
from os.path import isfile, join
import subprocess
import pandas as pd
from glob import glob
from Bio import Entrez, SeqIO


#%% define a function to make a list of target genbank files names
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
    gl9=pd.DataFrame({'strain':gl8,'NCBI ID':gl8, 'Pathotype':range(len(onlyfiles)-1)})
    path = directory2+'/StrainInformation.xlsx'
    StrainInformation = gl9.to_excel(path)
                

#%% retrives genome file names from a folder and returns an excell file with target genomes information
directory1 =r'/Users/omidard/Desktop/multispecies/Supplementary/genomes'
directory2 =r'/Users/omidard/Desktop/multispecies/Supplementary'
get_file_names (directory1,directory2)

#%%
# Load the information on the five strains we will be working with in this tutorial
StrainsOfInterest=pd.read_excel(directory2+'/StrainInformation.xlsx')
print(StrainsOfInterest)
#%%
#The Reference Genome is as Described in the Base Reconstruction; in these tutorials iML1515
referenceStrainID='reactome'
targetStrainIDs=list(StrainsOfInterest['NCBI ID'])

#%%
# define a function to gather information of the downloaded strains from the GenBank files
def get_strain_info(directory1):
    files = glob('%s/*.gbk'%directory1)
    strain_info = []
    
    for file in files:
        handle = open(file)
        record = SeqIO.read(handle, "genbank")
        for f in record.features:
            if f.type=='source':
                info = {}
                info['file'] = file
                info['id'] = file.split('\\')[-1].split('.')[0]
                for q in f.qualifiers.keys():
                    info[q] = '|'.join(f.qualifiers[q])
                strain_info.append(info)
    return pd.DataFrame(strain_info)
#%%
# information on the downloaded strain
get_strain_info(directory1)
print(pd.DataFrame(strain_info))
#%%
# define a function to parse the Genbank file to generate fasta files for both protein and nucleotide sequences
def parse_genome(id, type='prot', in_folder='genomes', out_folder='prots', overwrite=1):

    in_file = '%s/%s.gbk'%(in_folder, id)
    out_file='%s/%s.fa'%(out_folder, id)
    files =glob('%s/*.fa'%out_folder)
    
    if out_file in files and overwrite==0:
        print (out_file, 'already parsed')
        return
    else:
        print ('parsing %s'%id)
    
    handle = open(in_file)
    
    fout = open(out_file,'w')
    x = 0
    
    records = SeqIO.parse(handle, "genbank")
    for record in records:
        for f in record.features:
            if f.type=='CDS':
                seq=f.extract(record.seq)
                
                if type=='nucl':
                    seq=str(seq)
                else:
                    seq=str(seq.translate())
                    
                if 'locus_tag' in f.qualifiers.keys():
                    locus = f.qualifiers['locus_tag'][0]
                elif 'gene' in f.qualifiers.keys():
                    locus = f.qualifiers['gene'][0]
                else:
                    locus = 'gene_%i'%x
                    x+=1
                fout.write('>%s\n%s\n'%(locus, seq))
    fout.close()
#%%
# Generate fasta files for 5 strains of interest
directory3=r'/Users/omidard/Dsesktop/multispecies/Supplementary/prots'
directory4=r'/Users/omidard/Desktop/multispecies/Supplementary/nucl'
for strain in targetStrainIDs:
    parse_genome(strain, type='prot', in_folder=directory1, out_folder=directory3)
    parse_genome(strain, type='nucl', in_folder=directory1, out_folder=directory4)
#%%
#Also generate fasta files for the reference strain
parse_genome(referenceStrainID, type='nucl', in_folder=directory1, out_folder=directory4)
parse_genome(referenceStrainID, type='prots', in_folder=directory1, out_folder=directory3)
    
 #%%

def make_blast_db(id,folder=r'/Users/omidard/Desktop/multispecies/Supplementary/prots',db_type='prot'):
    import os
    
    out_file ='%s/%s.fa.pin'%(folder, id)
    files =glob('%s/*.fa.pin'%folder)
    
    if out_file in files:
        print (id, 'already has a blast db')
        return
    if db_type=='nucl':
        ext='fna'
    else:
        ext='fa'

    cmd_line='makeblastdb -in %s/%s.%s -dbtype %s' %(folder, id, ext, db_type)
    
    print ('making blast db with following command line...')
    print (cmd_line)
    os.system(cmd_line)

#%%
for strain in targetStrainIDs:
    make_blast_db(strain,folder=r'/Users/omidard/Desktop/multispecies/Supplementary/prots',db_type='prot')
make_blast_db(referenceStrainID,folder=r'/Users/omidard/Desktop/multispecies/Supplementary/prots',db_type='prot')
#%%
# define a function to run BLASTp
directory5=r'/Users/omidard/Desktop/multispecies/Supplementary/bbh'
def run_blastp(seq,db,in_folder=directory3, out_folder=directory5, out=None,outfmt=6,evalue=0.001,threads=1):
    import os
    if out==None:
        out='%s/%s_vs_%s.txt'%(out_folder, seq, db)
        print(out)
    
    files =glob('%s/*.txt'%out_folder)
    if out in files:
        print (seq, 'already blasted')
        return
    
    print ('blasting %s vs %s'%(seq, db))
    
    db = '%s/%s.fa'%(in_folder, db)
    seq = '%s/%s.fa'%(in_folder, seq)
    cmd_line='blastp -db %s -query %s -out %s -evalue %s -outfmt %s -num_threads %i' \
    %(db, seq, out, evalue, outfmt, threads)
    
    print ('running blastp with following command line...')
    print (cmd_line)
    os.system(cmd_line)
    return out
#%%
# define a function to get sequence length 

def get_gene_lens(query, in_folder=directory3):

    file = '%s/%s.fa'%(in_folder, query)
    handle = open(file)
    records = SeqIO.parse(handle, "fasta")
    out = []
    
    for record in records:
        out.append({'gene':record.name, 'gene_length':len(record.seq)})
    
    out = pd.DataFrame(out)
    return out
#%%
# define a function to get Bi-Directional BLASTp Best Hits
def get_bbh(query, subject, in_folder='bbh'):    
    
    #Utilize the defined protein BLAST function
    run_blastp(query, subject)
    run_blastp(subject, query)
    
    query_lengths = get_gene_lens(query, in_folder='prots')
    subject_lengths = get_gene_lens(subject, in_folder='prots')
    
    #Define the output file of this BLAST
    out_file = '%s/%s_vs_%s_parsed.csv'%(in_folder,query, subject)
    files=glob('%s/*_parsed.csv'%in_folder)
    
    #Combine the results of the protein BLAST into a dataframe
    print ('parsing BBHs for', query, subject)
    cols = ['gene', 'subject', 'PID', 'alnLength', 'mismatchCount', 'gapOpenCount', 'queryStart', 'queryEnd', 'subjectStart', 'subjectEnd', 'eVal', 'bitScore']
    bbh=pd.read_csv('%s/%s_vs_%s.txt'%(in_folder,query, subject), sep='\t', names=cols)
    bbh = pd.merge(bbh, query_lengths) 
    bbh['COV'] = bbh['alnLength']/bbh['gene_length']
    
    bbh2=pd.read_csv('%s/%s_vs_%s.txt'%(in_folder,subject, query), sep='\t', names=cols)
    bbh2 = pd.merge(bbh2, subject_lengths) 
    bbh2['COV'] = bbh2['alnLength']/bbh2['gene_length']
    out = pd.DataFrame()
    
    # Filter the genes based on coverage
    bbh = bbh[bbh.COV>=0.25]
    bbh2 = bbh2[bbh2.COV>=0.25]
    
    #Delineate the best hits from the BLAST
    for g in bbh.gene.unique():
        res = bbh[bbh.gene==g]
        if len(res)==0:
            continue
        best_hit = res.loc[res.PID.idxmax()]
        best_gene = best_hit.subject
        res2 = bbh2[bbh2.gene==best_gene]
        if len(res2)==0:
            continue
        best_hit2 = res2.loc[res2.PID.idxmax()]
        best_gene2 = best_hit2.subject
        if g==best_gene2:
            best_hit['BBH'] = '<=>'
        else:
            best_hit['BBH'] = '->'
        out=pd.concat([out, pd.DataFrame(best_hit).transpose()])
    
    #Save the final file to a designated CSV file
    out.to_csv(out_file)
#%%
# Execute the BLAST for each target strain against the reference strain, save results to 'bbh' i.e. "bidirectional best
# hits" folder to create
# homology matrix

for strain in targetStrainIDs:
    get_bbh(referenceStrainID,strain, in_folder='bbh')
#%%
#Load all the BLAST files between the reference strain and target strains

blast_files=glob('%s/*_parsed.csv'%'bbh')

for blast in blast_files:
    bbh=pd.read_csv(blast)
    print (blast,bbh.shape)
#%%
#Load the base reconstruction to designate the list of genes within the model
import cobra
import os
from os.path import join
data_dir = r'/Users/omidard/Desktop'
model = cobra.io.load_matlab_model(join(data_dir, "LBR3.mat"))
listGeneIDs=[]
for gene in model.genes:
    listGeneIDs.append(gene.id)
#%%
#Create 2 matrices of N, rows where N is the number of model genes and M columns where M is the number of target strains
#One matrix will be populated with the PID results from the blasts and another with the mapping of gene locus tags

ortho_matrix=pd.DataFrame(index=listGeneIDs,columns=targetStrainIDs)
geneIDs_matrix=pd.DataFrame(index=listGeneIDs,columns=targetStrainIDs)

#%%
#Parse through each blast file and acquire pertinent information for each matrix for each of the base reconstruction genes
for blast in blast_files:
    bbh=pd.read_csv(blast)
    listIDs=[]
    listPID=[]
    for r,row in ortho_matrix.iterrows():
        try:
            currentOrtholog=bbh[bbh['gene']==r].reset_index()
            listIDs.append(currentOrtholog.iloc[0]['subject'])
            listPID.append(currentOrtholog.iloc[0]['PID'])
        except:
            listIDs.append('None')
            listPID.append(0)
    for col in ortho_matrix.columns:
        if col in blast:
            ortho_matrix[col]=listPID
            geneIDs_matrix[col]=listIDs
#%%
# In this tutoriao, genes with a greater than 80% PID are considered present in the target strain genome 
# and consequently less than 80% are considered absent from the target strain genome
for column in ortho_matrix:
    ortho_matrix.loc[ortho_matrix[column]<=70.0,column]=0
    ortho_matrix.loc[ortho_matrix[column]>70.0,column]=1
#%%
#Define a function to generate FNA from the GBK files
def gbk2fasta(gbk_filename):
    faa_filename = '.'.join(gbk_filename.split('.')[:-1])+'.fna'
    input_handle  = open(gbk_filename, "r")
    output_handle = open(faa_filename, "w")

    for seq_record in SeqIO.parse(input_handle, "genbank") :
        print ("Converting GenBank record %s" % seq_record.id)
        output_handle.write(">%s %s\n%s\n" % (
               seq_record.id,
               seq_record.description,
               seq_record.seq))

    output_handle.close()
    input_handle.close()
#%%
#Define function to run the BLASTn
def run_blastn(seq, db,outfmt=6,evalue=0.001,threads=1):
    import os
    out = 'nucl/'+seq+'_vs_'+db+'.txt'
    seq = 'nucl/'+seq+'.fa'
    db = 'genomes/'+db+'.fna'
    
    cmd_line='blastn -db %s -query %s -out %s -evalue %s -outfmt %s -num_threads %i' \
    %(db, seq, out, evalue, outfmt, threads)
    
    print ('running blastn with following command line...')
    print (cmd_line)
    os.system(cmd_line)
    return out
#%%
# make nucleotide sequence databases 
for strain in targetStrainIDs:
    make_blast_db(strain,folder='genomes',db_type='nucl')

#%%
# convert genbank files to fna files for strains of interest
for strain in targetStrainIDs:
    gbk2fasta('genomes/'+strain+'.gbk')
#%%
# perform uni-directional BLASTn hit
genome_blast_res=[]
for strain in targetStrainIDs:
    res = run_blastn(referenceStrainID,strain)
    genome_blast_res.append(res)

#%%
#define a function to parse through the nucleotide BLAST results and form one matrix of all the results
def parse_nucl_blast(infile):
    cols = ['gene', 'subject', 'PID', 'alnLength', 'mismatchCount', 'gapOpenCount', 'queryStart', 'queryEnd', 'subjectStart', 'subjectEnd', 'eVal', 'bitScore']
    data = pd.read_csv(infile, sep='\t', names=cols)
    data = data[(data['PID']>70) & (data['alnLength']>0.7*data['queryEnd'])]
    data2=data.groupby('gene').first()
    return data2.reset_index()
#%%
# parse the nucleotide blast matrix 
na_matrix=pd.DataFrame()
for file in genome_blast_res:
    genes =parse_nucl_blast(file)
    name ='.'.join(file.split('_')[-1].split('.')[:-1])
    na_matrix = na_matrix.append(genes[['gene','subject','PID']])
na_matrix = pd.pivot_table(na_matrix, index='gene', columns='subject',values='PID')
#%%
na_matrix.head()

#%%
# define a function to extract the sequence from fna file 
def extract_seq(g, contig, start, end):
    from Bio import SeqIO
    handle = open(g)
    records = SeqIO.parse(handle, "fasta")
    
    for record in records:
        if record.name==contig:
            if end>start:
                section = record[start:end]
            else:
                section = record[end-1:start+1].reverse_complement()
                
            seq = str(section.seq)
    return seq
#%%
#Define updated matrices that will include genes based on sequence evidence that were missing due to lack of annotation
ortho_matrix_w_unannotated = ortho_matrix.copy()
geneIDs_matrix_w_unannotated = geneIDs_matrix.copy()

#%%
#Define matrix of the BLASTn results for all the pertinent model genes
nonModelGenes=[]
for g in na_matrix.index:
    if g not in listGeneIDs:
        nonModelGenes.append(g)

na_model_genes=na_matrix.drop(nonModelGenes)

#%%
#For each strain in the ortho_matrix, identify genes that meet threshold of SEQ similarity, but missing from
#annotated ORFS. Additionally, look at the sequence to ensure that these cases do not have early stop codons indicating
#nonfunctional even if the NA seqs meet the threshold

pseudogenes = {}

for c in ortho_matrix.columns:
    
    orfs = ortho_matrix[c]
    genes = na_model_genes[c]
    # All the Model Genes that met the BLASTp Requirements
    orfs2 = orfs[orfs==1].index.tolist()
    # All the Model Genes based off of BLASTn similarity above threshold of 80
    genes2 = genes[genes>=80].index.tolist()
    # By Definition find the genes that pass sequence threshold but were NOT in annotated ORFs:
    unannotated = set(genes2) -set(orfs2)
    
    # Obtain sequences of this list to check for premature stop codons:
    data = 'nucl/NC_000913.3_vs_%s.txt'%c
    cols = ['gene', 'subject', 'PID', 'alnLength', 'mismatchCount', 'gapOpenCount', 'queryStart', 'queryEnd', 'subjectStart', 'subjectEnd', 'eVal', 'bitScore']
    data = pd.read_csv(data, sep='\t', names=cols)
    #
    pseudogenes[c] = {}
    unannotated_data = data[data['gene'].isin(list(unannotated))]
    for i in unannotated_data.index:
        gene = data.loc[i,'gene']
        contig = data.loc[i,'subject'] 
        start = data.loc[i,'subjectStart']
        end = data.loc[i,'subjectEnd']
        seq = extract_seq('genomes/%s.fna'%c,contig, start, end)
        # check for early stop codons - these are likely nonfunctional and shouldn't be included
        if '*' in seq:
            print (seq)
            pseudogenes[c][gene]=seq
            # Remove the gene from list of unannotated genes
            unannotated-set([gene])
            
    
    print (c, unannotated)
    
    # For pertinent genes, retain those based off of nucleotide similarity within the orthology matrix and geneIDs matrix
    ortho_matrix_w_unannotated.loc[unannotated,c]=1
    for g in unannotated:
        geneIDs_matrix_w_unannotated.loc[g,c] = '%s_ortholog'%g
    

#%%
#Save the Presence/Absence Matrix and geneIDs Matrix for future use
ortho_matrix_w_unannotated.to_csv('ortho_matrix.csv')
geneIDs_matrix_w_unannotated.to_csv('geneIDs_matrix.csv')

#%%

#import package needed
import cobra
import pandas as pd
from cobra.io import load_json_model
from glob import glob
from cobra.manipulation.delete import delete_model_genes, remove_genes
import cobra
import os
from os.path import join

#%%
#import optlang
#optlang.glpk_interface.Configuration()
#Load the base E. coli reconstruction iML1515
#model = load_json_model('iML1515.json')
#model
model = cobra.io.load_matlab_model(os.path.join(r'/Users/omidard/Desktop', "LBR.mat"))
#%%
## Load the previously generated homology matrix for E. coli strains of interest
hom_matrix=pd.read_csv('ortho_matrix.csv')
hom_matrix=hom_matrix.set_index('Unnamed: 0')
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
    #s0001 is an artificial gene used in iML1515 for spontaneous reactions and as such has no homologs,
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
#load the geneID matrix from the notebook1 
models=glob('%s/*.json'%'Models')
geneIDs_matrix=pd.read_csv('geneIDs_matrix.csv')
geneIDs_matrix=geneIDs_matrix.set_index('Unnamed: 0')
geneIDs_matrix

#%%
#Utilize the geneIDs matrix to update the GPRs in each of the strain-specific models with the proper gene ID

from cobra.manipulation.modify import rename_genes

for mod in models:
    model=cobra.io.load_json_model(mod)
    for column in geneIDs_matrix.columns:
        if column in mod:
            currentStrain=column
    
    IDMapping=geneIDs_matrix[currentStrain].to_dict()
    IDMappingParsed = {k:v for k,v in IDMapping.items() if v != 'None'}
    
    rename_genes(model,IDMappingParsed)
    cobra.io.json.save_json_model(model,mod, pretty=False)

#%%

# gather the general information on the draft models
for strain in hom_matrix.columns:
    model=cobra.io.load_json_model(str('Models/'+strain+'.json'))
    print (model.id,'Number of Model Genes:',len(model.genes),'Number of Model Reactions:',len(model.reactions))























