#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 15:19:34 2022

@author: omidard
"""

"""
install and import packages
make directories
"""



import os
import pandas as pd
from glob import glob
from Bio import Entrez, SeqIO
import cobra
from os.path import join
from os.path import isfile, join
from os import listdir


#generate data directories
directories = ['reference_genome_dir', 'target_genome_dir', 'prots_dir', 'nucls_dir', 'bbh_dir' ,'present_absence_dir','initial_models_dir', 'output_models_dir', 'gapfilled_models_dir', 'blast_exe_dir','temp_files','ref_model_dir','nuc_blast_dir']
for i in directories:
    cmd_line = 'mkdir '+i
    os.system(cmd_line)


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


def dl_genome(id, folder='genomes'): # be sure get CORRECT ID
    files=glob('%s/*.gbk'%folder)
    out_file = '%s/%s.gbk'%(folder, id)

    if out_file in files:
        print (out_file, 'already downloaded')
        return
    else:
        print ('downloading %s from NCBI'%id)
        
    from Bio import Entrez
    Entrez.email = "omidard@biosustain.dtu.dk"     #Insert email here for NCBI
    handle = Entrez.efetch(db="nucleotide", id=id, rettype="gb", retmode="text")
    fout = open(out_file,'w')
    fout.write(handle.read())
    fout.close()


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


def make_blast_db(id,folder='prots',db_type='prot'):
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

    cmd_line='/usr/local/ncbi/blast/bin/makeblastdb -in %s/%s.%s -dbtype %s' %(folder, id, ext, db_type)
    
    print (' making blast db with following command line...')
    print (cmd_line)
    os.system(cmd_line)

def run_blastp(seq,db,in_folder='/Users/omidard/Desktop/multispecies/Supplementary/prots', out_folder='/Users/omidard/Desktop/multispecies/Supplementary/bbh', out=None,outfmt=6,evalue=0.001,threads=1):
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
    cmd_line='/usr/local/ncbi/blast/bin/blastp -db %s -query %s -out %s -evalue %s -outfmt %s -num_threads %i' %(db, seq, out, evalue, outfmt, threads)
    
    print ('running blastp with following command line...')
    print (cmd_line)
    os.system(cmd_line)
    return out


def get_gene_lens(query, in_folder='prots'):

    file = '%s/%s.fa'%(in_folder, query)
    handle = open(file)
    records = SeqIO.parse(handle, "fasta")
    out = []
    
    for record in records:
        out.append({'gene':record.name, 'gene_length':len(record.seq)})
    
    out = pd.DataFrame(out)
    return out


def get_bbh(query, subject, in_folder='bbh'):    
    
    #Utilize the defined protein BLAST function
    run_blastp(query, subject)
    run_blastp(subject, query)
    
    query_lengths = get_gene_lens(query, in_folder='/Users/omidard/Desktop/multispecies/Supplementary/prots')
    subject_lengths = get_gene_lens(subject, in_folder='/Users/omidard/Desktop/multispecies/Supplementary/prots')
    
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


def run_blastn(seq, db,outfmt=6,evalue=0.001,threads=1):
    import os
    out = '/Users/omidard/Desktop/multispecies/Supplementary/nucl/'+seq+'_vs_'+db+'.txt'
    seq = '/Users/omidard/Desktop/multispecies/Supplementary/nucl/'+seq+'.fa'
    db = '/Users/omidard/Desktop/multispecies/Supplementary/genomes/'+db+'.fna'
    
    cmd_line='/usr/local/ncbi/blast/bin/blastn -db %s -query %s -out %s -evalue %s -outfmt %s -num_threads %i' \
    %(db, seq, out, evalue, outfmt, threads)
    
    print ('running blastn with following command line...')
    print (cmd_line)
    os.system(cmd_line)
    return out



def parse_nucl_blast(infile):
    cols = ['gene', 'subject', 'PID', 'alnLength', 'mismatchCount', 'gapOpenCount', 'queryStart', 'queryEnd', 'subjectStart', 'subjectEnd', 'eVal', 'bitScore']
    data = pd.read_csv(infile, sep='\t', names=cols)
    data = data[(data['PID']>80) & (data['alnLength']>0.8*data['queryEnd'])]
    data2=data.groupby('gene').first()
    return data2.reset_index()


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


def cdm_fba(model,obj):
    parametr = 'BIOMASS'
    glc = 'EX_glc_D_e'
    arg = 'EX_arg_L_e'
    cys = 'EX_cys_L_e'
    glu = 'EX_glu_L_e'
    ile = 'EX_ile_L_e'
    leu = 'EX_leu_L_e'
    met = 'EX_met_L_e'
    phe = 'EX_phe_L_e'
    thr = 'EX_thr_L_e'
    tyr = 'EX_tyr_L_e'
    val = 'EX_val_L_e'
    cit = 'EX_cit_e'
    asc = 'EX_ascb_L_e'
    abz = 'EX_4abz_e'
    ade = 'EX_ade_e'
    fol = 'EX_fol_e'
    gly = 'EX_gly_e'
    gua = 'EX_gua_e'
    ins = 'EX_ins_e'
    ala = 'EX_ala_L_e'
    asp = 'EX_asp_L_e'
    his = 'EX_his_L_e'
    lys = 'EX_lys_L_e'
    pro = 'EX_pro_L_e'
    ser = 'EX_ser_L_e'
    trp = 'EX_trp_L_e'
    pyd = 'EX_pydam_e'
    pdx = 'EX_pydxn_e'
    rib = 'EX_ribflv_e'
    ac = 'EX_ac_e'
    thm = 'EX_thm_e'
    tym = 'EX_thymd_e'
    ura = 'EX_ura_e'
    cl = 'EX_cl_e'
    he = 'EX_h_e'
    h2o = 'EX_h2o_e'
    nh4 = 'EX_nh4_e'
    btn = 'EX_btn_e'
    ca2 = 'EX_ca2_e'
    k = 'EX_k_e'
    co = 'EX_co_e'
    co2 = 'EX_co2_e'
    pi = 'EX_pi_e'
    pent = 'EX_pnto_R_e'
    model.reactions.get_by_id(glc).lower_bound = -25
    model.reactions.get_by_id(arg).lower_bound = -10
    model.reactions.get_by_id(cys).lower_bound = -10
    model.reactions.get_by_id(glu).lower_bound = -10
    model.reactions.get_by_id(ile).lower_bound = -10
    model.reactions.get_by_id(leu).lower_bound = -10
    model.reactions.get_by_id(met).lower_bound = -10
    model.reactions.get_by_id(tyr).lower_bound = -10
    model.reactions.get_by_id(phe).lower_bound = -10
    model.reactions.get_by_id(thr).lower_bound = -10
    model.reactions.get_by_id(val).lower_bound = -10
    model.reactions.get_by_id(cit).lower_bound = -10
    model.reactions.get_by_id(asc).lower_bound = -10
    model.reactions.get_by_id(abz).lower_bound = -10
    model.reactions.get_by_id(ade).lower_bound = -10
    model.reactions.get_by_id(fol).lower_bound = -10
    model.reactions.get_by_id(gly).lower_bound = -10
    model.reactions.get_by_id(gua).lower_bound = -10
    model.reactions.get_by_id(ins).lower_bound = -10
    model.reactions.get_by_id(ala).lower_bound = -10
    model.reactions.get_by_id(asp).lower_bound = -10
    model.reactions.get_by_id(his).lower_bound = -10
    model.reactions.get_by_id(lys).lower_bound = -10
    model.reactions.get_by_id(pro).lower_bound = -10
    model.reactions.get_by_id(ser).lower_bound = -10
    model.reactions.get_by_id(trp).lower_bound = -10
    model.reactions.get_by_id(pyd).lower_bound = -1000
    model.reactions.get_by_id(pdx).lower_bound = -1000
    model.reactions.get_by_id(rib).lower_bound = -10
    model.reactions.get_by_id(ac).lower_bound = -10
    model.reactions.get_by_id(thm).lower_bound = -1000
    model.reactions.get_by_id(tym).lower_bound = -1000
    model.reactions.get_by_id(ura).lower_bound = -10
    model.reactions.get_by_id(cl).lower_bound = -1000
    model.reactions.get_by_id(he).lower_bound = -1000
    model.reactions.get_by_id(h2o).lower_bound = -1000
    model.reactions.get_by_id(nh4).lower_bound = -1000
    model.reactions.get_by_id(btn).lower_bound = -1000
    model.reactions.get_by_id(ca2).lower_bound = -1000
    model.reactions.get_by_id(k).lower_bound = -1000
    model.reactions.get_by_id(co).lower_bound = -1000
    model.reactions.get_by_id(co2).lower_bound = -1000
    model.reactions.get_by_id(pi).lower_bound = -1000
    model.reactions.get_by_id(parametr).upper_bound = 10
    model.reactions.get_by_id(parametr).lower_bound = -0
    model.reactions.get_by_id(pent).lower_bound = -1000
    model.objective = obj
    solution1 = model.optimize()
    print(solution1)
    model.summary()
    sol = str(solution1).split(None, 1)[1]
    solution = sol.split(None, 1)[0]
    return print(solution)



























