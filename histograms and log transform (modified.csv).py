#!/usr/bin/env python
# coding: utf-8

# In[12]:


import pandas as pd
from matplotlib import pyplot as plt
import numpy as np

#loading df
file = 'C:/Users/{kr.pA}/Downloads/methylated data/modified.csv'
df = pd.read_csv(file, sep = ',')

#editing df to desired format
df = df.transpose()

df.insert(416349, 416349, ['Patient Number','Patient 3 rep 1', 'Patient 3 rep 1','Patient 3 rep 1','Patient 3 rep 1','Patient 9','Patient 9','Patient 9','Patient 9','Patient 10','Patient 10',
             'Patient 11','Patient 11','Patient 11','Patient 11','Patient 12 rep 1','Patient 12 rep 1','Patient 14','Patient 14','Patient 14 rep 1','Patient 14 rep 1',
             'Patient 1','Patient 1','Patient 1','Patient 1','Patient 2','Patient 2','Patient 2','Patient 2','Patient 3','Patient 3',
             'Patient 3', 'Patient 3','Patient 4', 'Patient 4','Patient 4', 'Patient 4','Patient 7','Patient 7','Patient 7','Patient 7',
             'Patient 8','Patient 8','Patient 8','Patient 8','Patient 9','Patient 9','Patient 10','Patient 10','Patient 12','Patient 12',
             'Patient 14','Patient 14','Patient 15','Patient 15','Patient 15','Patient 15','Patient 16','Patient 16','Patient 16','Patient 16',
             'Patient 18','Patient 18','Patient 18','Patient 18','Patient 19','Patient 19','Patient 19','Patient 19','Patient 20','Patient 20',
             'Patient 20','Patient 20','Patient 21','Patient 21','Patient 21','Patient 21','Patient 22','Patient 22','Patient 22','Patient 22',
             'Patient 26','Patient 26','Patient 26','Patient 26','Patient 28','Patient 28','Patient 28','Patient 28','Patient 29','Patient 29',
             'Patient 29','Patient 29','Patient 29 rep 2','Patient 29 rep 2'], True) 

#populating a list of the genes and the patient number heading to rename the columns and initialize features list
genelist = list(df.iloc[0,0:])

features = genelist[:-1]
df.columns = genelist

#print(genelist)

df = df.drop(['ID_REF'], axis = 0)

#transform and create histograms out of nparray slicing

df_data = df.drop(['Patient Number'], axis =1)
npdata = df_data.to_numpy()

# need to separate out by methylated and unmethylated values
# we do this by deleting even and odd indices

meth = df_data.iloc[::2].to_numpy()
unmeth = df_data.iloc[1::2].to_numpy()


#separated histograms
print("untransformed methylated")
for i in meth:
    plt.hist(i, bins = 50)
plt.show()

print("untransformed unmethylated")
for i in unmeth:
    plt.hist(i, bins = 50)
plt.show()

#combined graphs
print("untransformed total")
for i in npdata:
    plt.hist(i, bins = 50)
plt.show()

print("log transformed total")
for i in npdata:
    plt.hist(i, bins = 50, log = True)
plt.show()
    

    

