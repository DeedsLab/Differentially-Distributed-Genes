import pandas as pd
import numpy as np
import os as os
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy.stats import variation
from scipy.stats import binom
import numpy.random as random
import scanpy as sc

large_root = "/media/timothyhamilton/data1/Tim_Hamilton/Limb_bud/"

# read the data from 10x .mtx:
print("Reading Data")
#skiprows=[i for i in range(11312,90637)
usec = list(range(0,11311))
newframe = pd.read_csv(large_root + "/limbbud.tsv", sep="\t", index_col = 0, usecols=usec)
metaframe = pd.read_csv(large_root + "/limbbud_meta_all.tsv", index_col = 0)
metaframe2 = metaframe.iloc[list(range(0,11310)), :]
cellheads = metaframe2.index.tolist()
newframe.columns = cellheads

clust_list = [newframe]
label = ['limbbud_13a']

for d, l in zip(clust_list, label):
    d.to_csv(large_root + "/" + l + ".csv")
    gene_headers = np.array(d.index)
    print("Cluster:" + l)
    data = d.to_numpy(dtype=float)
    avg_oall = np.mean(data, 1)
    sum_oall = np.sum(data, 1)
    sumoall = pd.DataFrame(sum_oall)
    sumoall.to_csv(large_root + "/" + l + "countsum_percell.csv", header=None, index=None)
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

