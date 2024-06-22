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


large_root = "/home/breannesparta/ddgs"
print("reading data")
label ="zheng10k"
data = pd.read_csv("/home/breannesparta/ddgs/10000_rep_1zheng9_resamp.csv", index_col = 0, nrows = 1)
#data = data.T
#it should be cell x gene, it is
print(data.head())

data1 = pd.read_csv(large_root + "/zheng_WC_r1_10k_genes_0_5kWilcox_test_results.csv", index_col = 0)
data2 = pd.read_csv(large_root + "/zheng_WC_r1_10k_genes_5k_9999Wilcox_test_results.csv", index_col = 0)
data3 = pd.read_csv(large_root + "/zheng_WC_r1_10k_genes_10k_14999Wilcox_test_results.csv", index_col = 0)
data4 = pd.read_csv(large_root + "/zheng_WC_r1_10k_genes_15k_19kWilcox_test_results.csv", index_col = 0)

#combine data into wilcox_frame for all genes
wilcox_frame = pd.concat([data1,data2,data3,data4], ignore_index=True) #this gives 4 rows with index 900; prob causing sort issue
#keeping original index in concat; causes but in sort
print(wilcox_frame.head())
print(wilcox_frame.shape)
wilcox_frame.to_csv(large_root + "/" + label + "Wilcox_test_results.csv")

cell_types = ['b','naive_cy','mono','regulatory','cd4','cd56','mem','naive_t','cytotoxic']
num_tests = wilcox_frame.shape[0]
#num_tests= len(gene_list)*len(cell_types) #is this the issue!???; WC frame is already multiplied
print(num_tests)

#BH correct after getting all gene statistics
print("Filtering N/As")
wilcox_frame=wilcox_frame.dropna()
print("Sorting and Finding Ranks")
#not sorting
print(wilcox_frame.loc[900,:])
sort_frame=wilcox_frame.sort_values(["Pvalue"])
print(sort_frame.loc[900,:])
sort_frame["Rank"]=sort_frame["Pvalue"].rank()
rank_list=sort_frame.iloc[:,4].values.tolist()
p_list=sort_frame.iloc[:,3].values.tolist()
adj_list = [(r*0.01)/num_tests for r in rank_list]
sort_frame.insert(5,"BH",adj_list)

## find largest pval that is smaller than iQm
sort_frame['issig'] = np.where(sort_frame['Pvalue'] < sort_frame['BH'], 1, 0)
sort_frame['lowv'] = np.where(sort_frame['issig'] == 1, sort_frame['Pvalue'], 0)
maxsig = np.max(sort_frame['lowv'])
print(maxsig)
sort_frame['sigp'] = np.where(sort_frame['Pvalue'] < maxsig, 1, 0)
#sgnames = wilcox_frame.index[wilcox_frame['sigp'] == 1].tolist()

#holder_list= list(set(subset_frame.iloc[:,0].values.tolist()))
#WCgenes = set(sgnames)
#sig_genes = pd.DataFrame(WCgenes)
#is the cut off value right?

max_index=0
cut_off_val=0
for i in range(len(adj_list)):
    if p_list[i]<adj_list[i] and i >= max_index:
        max_index=i
        cut_off_val= p_list[i]
print(cut_off_val)
print(max_index)
print(len(adj_list))

subset_frame=sort_frame.iloc[:max_index+1,:]
subset_frame.to_csv(large_root + "/" + label + "Wilcox_test_subset.csv") #has ranks and bh vals
holder_list= list(set(subset_frame.iloc[:,0].values.tolist()))
adj_data=data.loc[:,holder_list]
print("Printing")
adj_data.to_csv(large_root + "/" + label + "Wilcox.csv")
print("Done")