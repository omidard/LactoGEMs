#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 10 23:14:40 2022

@author: omidard
"""
#%%
import pandas as pd
import numpy as np
#%%
from Bio.KEGG import REST as kegg
from Bio.KEGG import Enzyme

def get_kegg_data(EC):
    kegg_output = kegg.kegg_get(EC).read()
    results = {}
    for line in kegg_output.split('\n'):
        splits = line.split()
        if not line.startswith(' '):    
            if len(splits) > 0:
                key = splits[0]
                value = ' '.join(splits[1:])
                results[key] = value
        else:
            results[key] += ' '.join(splits)
    return pd.DataFrame(results, index=[EC])


get_kegg_data_v = np.vectorize(get_kegg_data)

def get_kegg_info(kegg_ids):
    if isinstance(kegg_ids, str):
        kegg_ids = [kegg_ids]
    return pd.concat(get_kegg_data_v(kegg_ids), sort=False)
#%%
from Bio.KEGG import REST
import pandas as pd

request = get_kegg_data("ec:1.17.4.1")
keggsyn=list(request.NAME)
keggsyn2=str(keggsyn)
keggsyn3=keggsyn2.replace(";","\n")
keggsyn4=keggsyn3.replace('["','')
keggsyn5=keggsyn4.replace('"]','')