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

large_root = "/home/breanne/sourcedata/Zheng9"
print("reading data")
label ="zheng5k"
data= pd.read_csv(large_root + "/zheng5k_dropzeros.csv", index_col = 0)
data = data.T

data1 = pd.read_csv(large_root + "/genes1_5kWilcox_test_results.csv", index_col = 0)
data2 = pd.read_csv(large_root + "/genes5001_10kWilcox_test_results.csv", index_col = 0)
data3 = pd.read_csv(large_root + "/genes10001_15kWilcox_test_results.csv", index_col = 0)
data4 = pd.read_csv(large_root + "/genes15001_119324Wilcox_test_results.csv", index_col = 0)

#combine data into wilcox_frame for all genes
wilcox_frame = pd.concat([data1,data2,data3,data4])
wilcox_frame.to_csv(large_root + "/" + label + "Wilcox_test_results.csv")

cell_types = ['b','naive_cy','mono','regulatory','cd4','cd56','mem','naive_t','cytotoxic']
gene_list = wilcox_frame.keys()
num_tests= len(gene_list)*len(cell_types)

#BH correct after getting all gene statistics
print("Filtering N/As")
wilcox_frame=wilcox_frame.dropna()
print("Sorting and Finding Ranks")
wilcox_frame=wilcox_frame.sort_values(by="Pvalue")
wilcox_frame["Rank"]=wilcox_frame["Pvalue"].rank()
rank_list=wilcox_frame.iloc[:,4].values.tolist()
p_list=wilcox_frame.iloc[:,3].values.tolist()
adj_list = [r*0.01/num_tests for r in rank_list]
wilcox_frame.insert(5,"BH",adj_list)
max_index=0
cut_off_val=0
for i in range(len(adj_list)):
    if p_list[i]<adj_list[i] and i >= max_index:
        max_index=i
        cut_off_val= p_list[i]
print(cut_off_val)
#co .999996 ~ > 2k genes
subset_frame=wilcox_frame.iloc[:max_index+1,:]
subset_frame.to_csv(large_root + "/" + label + "Wilcox_test_subset.csv") #has ranks and bh vals

holder_list= list(set(subset_frame.iloc[:,0].values.tolist()))

adj_data=data.loc[:,holder_list]
print("Printing")
adj_data.to_csv(large_root + "/" + label + "Wilcox.csv")
print("Done")






