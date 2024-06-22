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
data_path = "/home/breannesparta/ddgs/10000_rep_1zheng9_resamp.csv"

print("reading data")
label ="zheng_WC_r1_10k_genes_15k_19k"
#skipi = list(range(1,1000)) + list(range(2001,19919))
skipi = list(range(1,14999))
cols = list(range(15000,19919))
cols.insert(0, 0)
data = pd.read_csv(data_path, index_col = 0, usecols=cols)
#data = data.T
#data= data.sort_values(by=data.keys()[0]) #cell #why

cell_types = ['b','naive_cy','mono','regulatory','cd4','cd56','mem','naive_t','cytotoxic']
print("Sorting by Cell Type") #this is the bug it includes and sorts by barcode
cell_types.sort()

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
        print("Doing "+g+" for cell type: " + c + "Coord: "+str((i,j)))
        # get coordinates of each gene per cell type, and print values into list

        #cellgroup1 = data.filter(regex=c+'.*', axis=0)
        cellgroup1=data[data.index.str.startswith(c)]

        first_list=cellgroup1.loc[:,g].values.tolist() # gene is in cell group c

        cellgroup2 = data.drop(cellgroup1.index)  # all cells except cell group c
        second_list = cellgroup2.loc[:,g].values.tolist() # gene values

        #first_list=data.loc[data[data.keys()[0]]== c].loc[:,g].values.tolist() # gene is in cell group c
        #second_list=data.loc[data[data.keys()[0]]!= c].loc[:,g].values.tolist() #gene is not in cell group c
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
wilcox_frame.to_csv(large_root + "/" + label + "Wilcox_test_results.csv")
"""
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
subset_frame=wilcox_frame.iloc[:max_index+1,:]
subset_frame.to_csv(large_root + "/" + label + "Wilcox_test_subset.csv")

holder_list= list(set(subset_frame.iloc[:,0].values.tolist()))


adj_data=data.loc[:,holder_list]
print("Printing")
adj_data.to_csv(large_root + "/" + label+ "Wilcox.csv")
print("Done")"""






