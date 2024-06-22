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



# read the data from 10x .mtx:
def neighbors(data, k=20):
    # for a given dataset, finds the k nearest neighbors for each point
    nbrs = NearestNeighbors(n_neighbors=k+1, algorithm='ball_tree').fit(data)
    distances, indices = nbrs.kneighbors(data)
    return indices[:,1:]

def jaccard(A,B):
    # for two sets A and B, finds the Jaccard distance J between A and B
    A = set(A)
    B = set(B)
    union = list(A|B)
    intersection = list(A & B)
    J = ((len(union) - len(intersection))/(len(union)))
    return(J)
def equalizer(A):
    for i in range(A.shape[1]):
        if np.mean(A[:,i])>=1.0:
            A[:,i]/=np.mean(A[:,i])
    return A



data_path = "/media/timothyhamilton/data1/Tim_Hamilton/Zheng_DDG/Non_Confusion/Full.csv"
#pheno_path = "/media/timothyhamilton/data1/Tim_Hamilton/R01_Ref_Datasets/Duo_"+d+"mix_id.csv"
print("reading data")
data = pd.read_csv(data_path)
data= data.sort_values(by=data.keys()[0])
data=data.set_index(data.keys()[0])
#data=data.transpose()
index_list=data.index.values.tolist()
cell_name_list = [i[:-16] for i in index_list]
#data.insert(0,"Cell Type",cell_name_list)
print("Sorting by Cell Type")

cell_types =list(set(cell_name_list))
cell_types.sort()

#cell_ids = data.index.values.tolist()


#cell_values = data.iloc[:,0].values.tolist()
print("Setting env")
gene_list = data.keys()
gene_col= [-1]*(len(gene_list)*len(cell_types))
type_col= [-1]*(len(gene_list)*len(cell_types))
stat_col= [-1]*(len(gene_list)*len(cell_types))
p_val_col= [-1]*(len(gene_list)*len(cell_types))

num_tests= len(gene_list)*len(cell_types)
indexer=0
for i,g in enumerate(gene_list):
    for j,c in enumerate(cell_types):
        print("Doing "+g+" for cell type: "+ c+ "Coord: "+str((i,j)))
        first_list=data.loc[data[data.keys()[0]]== c].loc[:,g].values.tolist()
        second_list=data.loc[data[data.keys()[0]]!= c].loc[:,g].values.tolist()
        try:
            stat,pval=ss.mannwhitneyu(first_list,second_list,alternative='two-sided')
        except ValueError:
            print("Possible Error Detected")
            stat=np.NaN
            pval=np.NaN
        stat_col[indexer]=stat
        gene_col[indexer]=g
        type_col[indexer]=c
        p_val_col[indexer]=pval
        indexer+=1

print("Doing Correction")
wilcox_frame=pd.DataFrame({"Gene":gene_col,"Cell_type":type_col,"Statistic":stat_col,"Pvalue": p_val_col})
wilcox_frame.to_csv("/media/timothyhamilton/data1/Tim_Hamilton/Zheng_DDG/Non_Confusion/Wilcox_test_results.csv")
print("Filtering N/As")
wilcox_frame=wilcox_frame.dropna()
print("Sorting and Finding Ranks")
wilcox_frame=wilcox_frame.sort_values(by="Pvalue")
wilcox_frame["Rank"]=wilcox_frame["Pvalue"].rank()
rank_list=wilcox_frame.iloc[:,4].values.tolist()
p_list=wilcox_frame.iloc[:,3].values.tolist()
adj_list = [r*0.05/num_tests for r in rank_list]
wilcox_frame.insert(5,"BH",adj_list)
max_index=0
cut_off_val=0
for i in range(len(adj_list)):
    if p_list[i]<adj_list[i] and i >= max_index:
        max_index=i
        cut_off_val= p_list[i]
print(cut_off_val)
subset_frame=wilcox_frame.iloc[:max_index+1,:]
subset_frame.to_csv("/media/timothyhamilton/data1/Tim_Hamilton/Zheng_DDG/Non_Confusion/Wilcox_test_results_subset.csv")
holder_list= list(set(subset_frame.iloc[:,0].values.tolist()))


#collection_array=np.asarray(collection_array)
adj_data=data.loc[:,holder_list]
print("Printing")
adj_data.to_csv("/media/timothyhamilton/data1/Tim_Hamilton/Zheng_DDG/Non_Confusion/Wilcox.csv")
print("Done")
#cell_ids_frame= pd.DataFrame({"Barcode": cell_ids})
#gene_frame= pd.DataFrame(index=gene_list)





