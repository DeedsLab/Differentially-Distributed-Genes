import pandas as pd
import numpy as np
import os as os
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy.stats import variation
from scipy.stats import binom
import numpy.random as random
import scanpy as sc

large_root = "/media/timothyhamilton/data1/Tim_Hamilton/Limb_bud"

# read the data from 10x .mtx:
print("Reading Data")
#skiprows=[i for i in range(11312,90637)
usec = list(range(0,11311))
newframe = pd.read_csv(large_root + "/limbbud.csv", index_col = 0, usecols=usec)
metaframe = pd.read_csv(large_root + "/limbbud_meta_all.csv", index_col = 0)
metaframe2 = metaframe.iloc[list(range(0,11310)), :]
cellheads = metaframe2.index.tolist()
newframe.columns = cellheads

##now we un normalize
#x = ncounts / 10000
#count = x*normalized value

countfactor = [x / 10000 for x in metaframe2.n_counts]
newframe = newframe.T
#i = countfactor.shape[0]
#j = newframe.shape[0]
#print("metadata"+str(i))
#print("count_cell_data"+str(j))
newframe.mul(countfactor, axis =0) #error length must be 43346 given 11310 where does 43 come from
newframe = newframe.T

dropzero = newframe[np.all(newframe == 0, axis=1)].index
newframe = newframe.drop(dropzero, 'index') #drop genes that are not detected in any cells
#newframe.to_csv(large_root + "/limbbud_12_13.csv")

#cellgroup1 = newframe.filter(regex='^limb8_15',axis=1) #b_cells
#cellgroup2 = newframe.filter(regex='^limb7_10',axis=1) #cytotoxic
#cellgroup3 = newframe.filter(regex='^limb6_15',axis=1) #monocytes
#cellgroup4 = newframe.filter(regex='^limb5_13',axis=1) #t
#cellgroup5 = newframe.filter(regex='^limb4_12',axis=1) #thelper
#cellgroup6 = newframe.filter(regex='^limb3_11',axis=1) #nk
#cellgroup7 = newframe.filter(regex='^limb1_13',axis=1) #memory t
#cellgroup8 = newframe.filter(regex='^limb13_14',axis=1) #naive t
#cellgroup9 = newframe.filter(regex='^limb12_13',axis=1) #cytotoxic_t

#clust_list = [cellgroup1, cellgroup2, cellgroup3, cellgroup4, cellgroup5,cellgroup6, cellgroup7, cellgroup8, cellgroup9]
#label = ['limbbud_15b','limbbud_10_5','limbbud_15a','limbbud_13b','limbbud_12','limbbud_11','limbbud_13_5','limbbud_14','limbbud_13a']

clust_list = [newframe]
label = ['limbbud_13a_denorm']

for d, l in zip(clust_list, label):
    d.to_csv(large_root + "/" + l + ".csv")
    gene_headers = np.array(d.index)
    print("Cluster:" + l)
    data = d.to_numpy(dtype=float)
    avg_oall = np.mean(data, 1)
    data[data == 0] = np.nan
    num_cells = np.count_nonzero(~np.isnan(data),1) #how many cells each mRNA is detected in
    NT = data.shape[1] #total ncells there are 17k genes in hydra! 25k cells

    numc = str(NT)
    getp = np.vstack((gene_headers, num_cells, avg_oall))
    getpar = np.transpose(getp)
    getparams = pd.DataFrame(getpar)
    getparams.to_csv(large_root + "/" + l + "gene_params.csv", header=None, index=None)
    base_fn = (l + numc + "gene_params.txt")
    with open(os.path.join(large_root, base_fn),'w') as outfile:
        getparams.to_string(outfile, header=None, index=None)

