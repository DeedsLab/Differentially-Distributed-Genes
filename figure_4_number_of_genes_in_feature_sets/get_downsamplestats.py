import pandas as pd
import numpy as np
import os as os
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy.stats import variation
from scipy.stats import binom
import numpy.random as random
import scanpy as sc

large_root = "/home/breannesparta/ddgs/zheng9/cellgene"
rep_list = [1,2,3,4,5,6,7,8,9,10]
# read the data from 10x .mtx:
print("Reading Data ")

#cols = [c for c in new_frame.columns if c[:4] != 'Cell']
#new_frame = new_frame[cols]
filename_list = ['zheng5kdownsample0.9reps1','zheng5kdownsample0.9reps2','zheng5kdownsample0.9reps3','zheng5kdownsample0.9reps4','zheng5kdownsample0.9reps5','zheng5kdownsample0.9reps6','zheng5kdownsample0.9reps7','zheng5kdownsample0.9reps8','zheng5kdownsample0.9reps9','zheng5kdownsample0.9reps10']
label = ['z5k0_09ds_rep1_','z5k0_09ds_rep2_','z5k0_09ds_rep3_','z5k0_09ds_rep4_','z5k0_09ds_rep5_','z5k0_09ds_rep6_','z5k0_09ds_rep7_','z5k0_09ds_rep8_','z5k0_09ds_rep9_','z5k0_09ds_rep10_']

for filename, l in zip(filename_list, label):
    new_frame = pd.read_csv(large_root + "/" + filename + ".csv", index_col=0) #cell by gene
    new_frame = new_frame.T
    dropzero = new_frame[np.all(new_frame == 0, axis=1)].index #gene x cell
    newframe = new_frame.drop(dropzero, 'index')
    d = new_frame

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
    #getparams.to_csv(large_root + "/" + l + "gene_params.csv", header=None, index=None)
    base_fn = (l + numc + "gene_params.txt")
    with open(os.path.join(large_root, base_fn),'w') as outfile:
        getparams.to_string(outfile, header=None, index=None)

