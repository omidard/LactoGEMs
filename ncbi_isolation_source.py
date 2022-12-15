#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  4 17:57:53 2022

@author: omidard
"""

#collect isolation source from NCBI
#requires informations from all_models_metabolite_production.py





#imports
from Bio import Entrez
import pandas as pd





isol=[]
gid=[]

for i in info.columns:
    target = i.replace('.json','')
    acc = target
    gid.append(target)
    search = Entrez.read(Entrez.esearch(db="assembly", term=acc))
    summary = Entrez.read(Entrez.esummary(db="assembly",id=search['IdList'],report="full"))
    bios = str(summary.items())
    biosa = bios.split("Biosource")[0]
    biosam = biosa.split("BioSampleId': '")[1]
    biosamp = biosam.replace("', '","")
    biosampl = biosamp.replace('"','')
    handle=Entrez.esummary(db='biosample',id=biosampl)
    recs=Entrez.read(handle)
    records = recs['DocumentSummarySet']['DocumentSummary'][0]
    if 'isolation source' in records['SampleData']:
        isolation_source = records['SampleData'].split('display_name="isolation source">')[1].split('</Attribute>')[0]
        isol.append(isolation_source)
    if 'isolation source' not in records['SampleData']:
        isolation_source='Not reported'
        isol.append(isolation_source)
    
isol_data = pd.DataFrame({'Genome_id':gid,'Isolation_source':isol})
isol_data