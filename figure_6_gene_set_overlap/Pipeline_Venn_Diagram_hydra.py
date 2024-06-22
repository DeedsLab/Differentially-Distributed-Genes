import pandas as pd
import numpy as np
import numpy.random as random
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import anndata as ad
from sklearn.neighbors import NearestNeighbors
from sklearn.neighbors import DistanceMetric
from venn import venn
from matplotlib_venn import venn2, venn3

import smtplib
import os
import sys

#HVG_frame= pd.read_csv("/media/timothyhamilton/data1/Tim_Hamilton/Zheng_DDG/Non_Confusion/HVGs.csv",nrows=1)
#HVG_set= set(HVG_frame.keys()[1:])
DDG_frame=pd.read_csv("/media/timothyhamilton/data1/Tim_Hamilton/Hydra/Hydra_DDG_data.csv",index_col = 0, usecols = [0,1])
DDG_set= set(DDG_frame.index)

Standard_HVG_Frame= pd.read_csv("/media/timothyhamilton/data1/Tim_Hamilton/Hydra/Hydra_HVGs.csv",nrows=1)
all_genes = Standard_HVG_Frame.keys()[1:]
newlist = []
for gene in all_genes:
    newgene = str(gene)
    a = newgene.replace('.','').replace('-','')
    newlist.append(a)
all_headers = np.array(newlist)
Standard_HVG_set= set(all_headers)

Drop_Frame= pd.read_csv("/media/timothyhamilton/data1/Tim_Hamilton/Hydra/Hydra_M3Drop_Genes.csv",nrows=1)
M3Drop_set= set(Drop_Frame.keys()[1:])

Townes_Frame= pd.read_csv("/media/timothyhamilton/data1/Tim_Hamilton/Hydra/Hydra_Raw_Cell_Types_Townes.csv",nrows=1)
Townes_set= set(Townes_Frame.keys()[1:])

title_list = [ "DDGs","HVGs","M3Drop","Townes"]
val_list = [DDG_set, Standard_HVG_set,M3Drop_set, Townes_set]
val_dict= dict(list(zip(title_list,val_list)))
fig = plt.figure(figsize=(15, 8))
ax1 = fig.add_subplot(111)
ax1.set_title("Overlap of Gene_Sets")
venn(val_dict, cmap="viridis", fontsize=12, legend_loc="upper left", ax=ax1)
ax1.set_rasterized(True)
plt.savefig("/media/timothyhamilton/data1/Tim_Hamilton/Hydra/Hydra_gene_overlap.eps")
plt.savefig("/media/timothyhamilton/data1/Tim_Hamilton/Hydra/Hydra_gene_overlap.png")
