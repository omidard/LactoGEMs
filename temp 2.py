# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
#import required mudules
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

#dict file to excell (!!!!because of size limiting it isnt used here!!!!)
mport pandas as pd
df = pd.DataFrame(data=cl, index=[0])
df = (df.T)
print (df)
df.to_excel('dict1.xlsx')


#load excel file
import pandas as pd
gpr = pd.read_excel ('ugpr.xls')
print (gpr)

#turn ugpr dataframe to a list
gprlist = gpr.values.tolist()

#%%
#transpose column to row (!!!!isnt used here!!!!)
gpr2_tr = gpr2.transpose()

#%%
#turn a list into a single string
listToStr = ' '.join([str(elem) for elem in gpr2.GPR])


#TURN STINGS WITH SPACE INTO A LIST
gpr3 = listToStr.split()

#%%
#Find common elements of set and list
list1_as_set = set(gpr3)
intersection = list1_as_set.intersection(genmap)
intersection_as_list = list(intersection)

#get values(cdhits) for a list of keys(locustags)

gennames = [genmap[x] for x in intersection_as_list]

#%%get keys(locustags) for values(cdhits)


def getKeysByValue(dictOfElements, valueToFind):
    listOfKeys = list()
    listOfItems = dictOfElements.items()
    for item  in listOfItems:
        if item[1] == valueToFind:
            listOfKeys.append(item[0])
    return  listOfKeys
for i in gennames:      ------------ should be modified
listOfKeys = getKeysByValue(genmap, valuetofind)
for key  in listOfKeys:
        print(key)


#%%
#Find the indices at which any element of one list occurs in another
st = set(intersection_as_list)
indices_list = [i for i, e in enumerate(cl) if e in st]


#extract cdhits using mapped indices
cdhits = list(genmap.items())


#extract items using indices list

gholi= [None] * len(indices_list)
m=0
for i in indices_list:
   a = cdhits[i]
   gholi[m]=a
   m= m+1
   print(m)
  
#dict values to list
lclvl = list(genmap.values())



