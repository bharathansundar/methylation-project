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

from sklearn.preprocessing import StandardScaler

# Separating out the features
x = df.loc[:, features].values
#print(x)

# Separating out the target (patient # vals)
y = df.loc[:,['Patient Number']].values
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
df_Patient_Number = df[['Patient Number']]

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






