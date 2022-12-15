import cobra
import pandas as pd
from cobra.io import load_json_model
from glob import glob
from cobra.manipulation.delete import delete_model_genes, remove_genes
import os
from os.path import join
import pandas as pd
from glob import glob
from Bio import Entrez, SeqIO
from os.path import join
from os.path import isfile, join
from os import listdir
#%%

flux2=[]
flux3 = []
for i in range(len(fluxes)):
    if fluxes[i] > 0:
        fx = fluxes.index[i]
        fx2 = fluxes[i]
        flux2.append(fx)
        flux3.append(fx2)
    if fluxes[i] < 0:
        fx = fluxes.index[i]
        fx2 = fluxes[i]
        flux2.append(fx)
        flux3.append(fx2)
fluxes2 = pd.DataFrame({'id':flux2,'flux':flux3})
#%%
for i in range(len(fluxes2.id)):
    if 'ATPS3r' in fluxes2.id[i]:
        print(fluxes2.flux[i])

#%%
unireactions = []
universal = cobra.io.load_json_model(join('gapfilled_models_dir', "model.json"))
for reaction in universal.reactions:
    tar = reaction.id
    unireactions.append(tar)


#%%
ov =[]
models=glob('%s/*.json'%'output_models_dir')
for model in models:
    model=cobra.io.load_json_model(model)
    mre=[]
    #col = model.id
    for reaction in model.reactions:
        que = reaction.id
        mre.append(que)
    ov.append(mre)
#%%
overal = pd.DataFrame()
for x in unireactions:
    counter = []
    for i in range(len(ov)):
        if x in ov[i]:
            counter.append(1)
        else:
            counter.append(0)
    overal[x] = counter
#%%
columns = list(overal)
#%%
ctc = []
for x in columns:
    ct = []
    for i in overal[x]:
        if i == 1:
            ct.append(1)
    ctc.append(ct)
#%%
fr = []
for i in range(len(ctc)):
    freq = int(len(ctc[i]))/133 * 100
    fr.append(freq)
#%%
leng =pd.DataFrame()
leng['id'] = columns
leng['freq'] = fr
#%%
leng2 = leng.sort_values('freq')

            
            
            
            
    
            
        

        