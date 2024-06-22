import pandas as pd
import numpy as np
import numpy.random as random
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import anndata as ad
import scipy.stats as ss
from sklearn.neighbors import NearestNeighbors
from sklearn.neighbors import DistanceMetric

import smtplib
import os
import sys


data_path = "/home/breanne/sourcedata/Zheng9/zheng5k_dropzeros.csv"
#pheno_path = "/media/timothyhamilton/data1/Tim_Hamilton/R01_Ref_Datasets/Duo_"+d+"mix_id.csv"
print("reading data")
data = pd.read_csv(data_path, index_col = 0)
data = data.T
index_list=data.index.values.tolist()
cell_name_list = [i[:-16] for i in index_list]
data.insert(0,"Cell Type",cell_name_list)
print("Sorting by Cell Type")
data= data.sort_values(by=data.keys()[0])
cell_types =list(set(cell_name_list))
cell_types.sort()

#cell_ids = data.index.values.tolist()
#cell_values = data.iloc[:,0].values.tolist()
print("Setting env")
gene_list = data.keys()[1:]
collection_array=[[-1]*len(cell_types) for i in range(len(gene_list))]
stat_array = np.empty((len(gene_list),len(cell_types)))
pval_array = np.empty((len(gene_list),len(cell_types)))
#collection_array[:]=np.NaN
stat_array[:]=0
pval_array[:]=0
cutoff= 0.01/(len(gene_list)*len(cell_types))
for i,g in enumerate(gene_list):
    for j,c in enumerate(cell_types):
        print("Doing "+g+" for cell type: "+ c+ "Coord: "+str((i,j)))
        first_list=data.loc[data[data.keys()[0]]== c].loc[:,g].values.tolist()
        second_list=data.loc[data[data.keys()[0]]!= c].loc[:,g].values.tolist()
        try:
            stat,pval=ss.mannwhitneyu(first_list,second_list,alternative='two-sided')
        except ValueError:
            stat=np.NaN
            pval=np.NaN
        stat_array[i,j]=stat
        pval_array[i,j]=pval
        if pval<=cutoff:
            collection_array[i][j]=g
print("Formatting for Printing")
holder_list=[]
for i in range(len(collection_array)):
    for j in range(len(collection_array[i])):
        label= collection_array[i][j]
        if label != -1 and label not in holder_list:
            holder_list.append(label)
#holder_list=list(set(holder_list))
#collection_array=np.asarray(collection_array)
adj_data=data.loc[:,holder_list]
print("Printing")
adj_data.to_csv("/home/breanne/sourcedata/Zheng9/zheng5k/Wilcox.csv")
pd.DataFrame(collection_array,index =gene_list,columns=cell_types).to_csv("/home/breanne/sourcedata/Zheng9/zheng5k/Wilcox_gene_list.csv")
pd.DataFrame(stat_array,index=gene_list,columns=cell_types).to_csv("/home/breanne/sourcedata/Zheng9/zheng5k/Wilcox_stat_array.csv")
pd.DataFrame(pval_array,index=gene_list,columns=cell_types).to_csv("/home/breanne/sourcedata/Zheng9/zheng5k/Wilcox_pval_array.csv")
print("Done")
