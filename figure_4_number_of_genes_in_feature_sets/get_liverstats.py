import pandas as pd
import numpy as np
import os as os
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy.stats import variation
from scipy.stats import binom
import numpy.random as random
import scanpy as sc

large_root = "/home/breanne/gcc_data/parse_liver"

# read the data from 10x .mtx:
print("Reading Data ")
new_frame = pd.read_csv(large_root + "/parse_liver.csv", index_col = 0)
cols = [c for c in new_frame.columns if c[:4] != 'Cell']
new_frame = new_frame[cols]
new_frame = new_frame.T
dropzero = new_frame[np.all(new_frame == 0, axis=1)].index
newframe = new_frame.drop(dropzero, 'index')
clust_list = [newframe]
label = ['parse_liver']

for d, l in zip(clust_list, label):
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

