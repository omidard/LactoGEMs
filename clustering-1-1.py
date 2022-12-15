#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# clustring genes related to each reaction (by gene name) - approved!

import pandas as pd
from Bio import SeqIO
from Bio import AlignIO
from Bio import Phylo
import os
import re

gpr = pd.read_excel (r'/Users/omidard/Desktop/reactome/gpr.xls') #xls list of gprs and related rxn_ids
fasta_file = r'/Users/omidard/Desktop/reactome/pan.faa'  # Input fasta file

gprdf = pd.DataFrame(gpr)
re_ge = dict(tuple(gprdf.groupby('id')))
r_g = list(re_ge.values())

fasta_file = r'/Users/omidard/Desktop/reactome/pan.faa'
filenamegenerator = list(re_ge.keys())
joker = list(range(len (re_ge)))
directoryin = r'/Users/omidard/Desktop/reactome/groupedgprsgenes/'
directoryout = r'/Users/omidard/Desktop/reactome/reactionsclusters/'
for i in joker:
    clusname = filenamegenerator[i]
    wanted_file = directoryin+'targetseqs'+str(filenamegenerator[i])+'.fasta'
    tehran = list(SeqIO.parse(wanted_file, "fasta"))
    ker=[]
    for i in range(len(tehran)):
        cph=tehran[i].description
        ker.append(cph)
        ods = str(ker)
        bud = ods.replace("[NADPH]","NADPH")
        lux = bud.replace("[glutamine-hydrolyzing]","glutamine-hydrolyzing")
        lux1 = lux.replace("[NADP(+)]", "NADP")
        lux2 = lux1.replace("[NAD(P)+]", "NADP-NAD")
        lux3 = lux2.replace("[acyl-carrier-protein]","ACP")
        lux4 = lux3.replace("(5-phosphoribosyl)-5-[(5-phosphoribosylamino)methylideneamino]", "5-phosphoribosyl-5-5-phosphoribosylamino methylideneamino")
        lux5 = lux4.replace("[NADH]", "NADH")
        man = lux5.replace("['","")
        cal = man.replace("']","")
        was = cal.replace("', '", "\n")
        gaz = re.findall(r'\[.*?\]', was)
        osl = str(gaz)
        maz = osl.replace("', '","\n")
        shi = maz.replace("['","")
        gil = shi.replace("']","") #final gene names
        san = re.sub("[\[].*?[\]]", "", was)
        jav = san.replace("  ", "\n") #final locoustags name
        for line in jav:
            bag = list(jav.splitlines()) #locustag list
        for line in gil:
            mos = list(gil.splitlines()) #gene name list
        earth = {'Locustags':bag, 'Gene_name':mos}
        mway = pd.DataFrame(earth)
        amd = list(tuple(mway.groupby('Gene_name')))
        for i in range(len(amd)):
            gmlt = amd[i][1]['Locustags']
            result_file = directoryout+'cluster'+str(clusname)+'_'+str(i)+'.fasta'
            with open(result_file, "w") as mde:
                for seq in tehran:
                    if seq.id in str(gmlt):
                        SeqIO.write([seq], mde, "fasta")
                        
                        
#%%
merge-gbk-records 1.gbk 2.gbk 3.gbk 4.gbk 5.gbk 6.gbk 7.gbk 8.gbk 9.gbk 10.gbk 11.gbk 12.gbk 13.gbk 14.gbk 15.gbk 16.gbk 17.gbk 18.gbk 19.gbk 20.gbk 21.gbk 22.gbk 23.gbk 24.gbk 25.gbk 26.gbk 27.gbk 28.gbk 29.gbk 30.gbk 31.gbk 32.gbk 33.gbk 34.gbk 35.gbk 36.gbk 37.gbk 38.gbk 39.gbk 40.gbk 41.gbk 42.gbk 43.gbk 44.gbk 45.gbk 46.gbk 47.gbk 48.gbk 49.gbk > merged.gbk






        