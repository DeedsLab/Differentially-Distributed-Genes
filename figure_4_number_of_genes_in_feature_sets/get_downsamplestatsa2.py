import pandas as pd
import numpy as np
import os as os
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy.stats import variation
from scipy.stats import binom
import numpy.random as random
import scanpy as sc

large_root = "/home/breannesparta/ddgs/A20_3T3/cellgene"
rep_list = [1,2,3,4,5,6,7,8,9,10]
# read the data from 10x .mtx:
print("Reading Data ")

#cols = [c for c in new_frame.columns if c[:4] != 'Cell']
#new_frame = new_frame[cols]

filename_list = ['a20_3t3downsample0.5reps1','a20_3t3downsample0.5reps2','a20_3t3downsample0.5reps3','a20_3t3downsample0.5reps4','a20_3t3downsample0.5reps5','a20_3t3downsample0.5reps6','a20_3t3downsample0.5reps7','a20_3t3downsample0.5reps8','a20_3t3downsample0.5reps9','a20_3t3downsample0.5reps10']
label = ['a203t3_05ds_rep1_','a203t3_05ds_rep2_','a203t3_05ds_rep3_','a203t3_05ds_rep4_','a203t3_05ds_rep5_','a203t3_05ds_rep6_','a203t3_05ds_rep7_','a203t3_05ds_rep8_','a203t3_05ds_rep9_','a203t3_05ds_rep10_']

for filename, l in zip(filename_list, label):
    new_frame = pd.read_csv(large_root + "/" + filename + ".csv", index_col=0) #cell by gene
    new_frame = new_frame.T
    dropzero = new_frame[np.all(new_frame == 0, axis=1)].index #gene x cell
    newframe = new_frame.drop(dropzero, 'index')
    d = new_frame
    print(d)

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

