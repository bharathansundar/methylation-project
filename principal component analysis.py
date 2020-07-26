#!/usr/bin/env python
# coding: utf-8

# In[18]:


import pandas as pd
from matplotlib import pyplot as plt

#loading df
file = 'C:/Users/{kr.pA}/Downloads/methylated data/modified.csv'
df = pd.read_csv(file, sep = ',')

#editing df to desired format
df = df.transpose()

#log transform
a = df.iloc[0:1]

b = df.iloc[1:]
b = b.astype(float)
b = b.transform(np.log)

df = pd.concat([a, b])

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

print(genelist)

df = df.drop(['ID_REF'], axis = 0)


#***SEPARATING PRIMARY VS RECURRENT***

missing = df.iloc[10:16]

#recur1 = unmeth
#recur2 = meth
#prim1 = unmeth
#prim2 = meth

recurmissing = missing.iloc[2:]

#recurrent

recur1 = df.iloc[3::4]
recur1 = recur1.drop("Patient 11_primary tumour Unmethylated signal")

recur1 = pd.concat([recur1, recurmissing.iloc[1::2]])

recur2 = df.iloc[2::4]
recur2 = recur2.drop("Patient 11_primary tumour Methylated signal")


s = df.iloc[-1:]
recur1 = pd.concat([recur1, s])

xs = df.iloc[-2:-1]
recur2= pd.concat([recur2, xs])
recur2= pd.concat([recur2, recurmissing.iloc[::2]])


#primary

prim1 = df.iloc[1::4]
prim1 = prim1.drop("Patient 11_recurrent tumour Unmethylated signal")
prim1 = prim1.drop("Patient 12_recurrent tumour Unmethylated signal")
prim1 = prim1.drop("Patient 29_recurrent tumour_2 Unmethylated signal")

prim1 = pd.concat([prim1, missing.iloc[1:2]])

prim2 = df.iloc[::4]
prim2 = prim2.drop("Patient 11_recurrent tumour Methylated signal")
prim2 = prim2.drop("Patient 12_recurrent tumour Methylated signal")
prim2 = prim2.drop("Patient 29_recurrent tumour_2 Methylated signal")

prim2= pd.concat([prim2,missing.iloc[0:1]])

#final indexed dataframes for the primary and recurrent data (separated by meth and unmeth within)

primary = pd.concat([prim1,prim2])   
recurrent = pd.concat([recur1, recur2])



from sklearn.preprocessing import StandardScaler

# Currently for transformed primary tumor data (can change to df, recurring, for full or recurrent data only)
# Separating out the features
x = primary.loc[:, features].values

# Separating out the target (patient # vals)
y = primary.loc[:,['Patient Number']].values
#print("Patient Number equals \n", y)

# Standardizing the features
x = StandardScaler().fit_transform(x)

#PCA starts here ***
from sklearn.decomposition import PCA

pca = PCA(n_components=2)
principalComponents = pca.fit_transform(x)
principalDf = pd.DataFrame(data = principalComponents
             , columns = ['principal component 1', 'principal component 2'])

# Extract the patient number column
df_Patient_Number = primary[['Patient Number']]

# Copy the index from principalDf to df_Patient_Number.
# If we don't do this, the pd.concat() will not know which row to match against which.
# So, it will end up doubling the final output, with nonsense (NaN) data all over the place.
df_Patient_Number.index = principalDf.index

finalDf = pd.concat([principalDf, df_Patient_Number], axis = 1)

# visualizing the data in 2d

fig = plt.figure(figsize = (300,300))
ax = fig.add_subplot(1,1,1) 
ax.set_xlabel('Principal Component 1', fontsize = 140)
ax.set_ylabel('Principal Component 2', fontsize = 140)
ax.set_title('2 component PCA', fontsize = 140)

patients = ['Patient 1','Patient 2','Patient 3','Patient 3 rep 1','Patient 4','Patient 7','Patient 9','Patient 10','Patient 11','Patient 12','Patient 12 rep 1','Patient 14','Patient 14 rep 1','Patient 15','Patient 16','Patient 18','Patient 19','Patient 20','Patient 21','Patient 22','Patient 26','Patient 28','Patient 29', 'Patient 29 rep 2']
colors = ['r', 'g','b','y','c','m','k','brown','gray','purple','teal','coral','goldenrod','pink','crimson','fuchsia','dodgerblue','chocolate','mediumpurple','olivedrab','chartreuse','darkslateblue','maroon','hotpink']


for patient, color in zip(patients,colors):
    indicesToKeep = finalDf['Patient Number'] == patient
    ax.scatter(finalDf.loc[indicesToKeep, 'principal component 1']
               , finalDf.loc[indicesToKeep, 'principal component 2']
               , c = color
               , s = 3000)
ax.legend(patients, fontsize = 110)
ax.grid()






