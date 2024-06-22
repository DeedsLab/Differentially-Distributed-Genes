import pandas as pd
import numpy as np
import os as os
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy.stats import variation
from scipy.stats import binom
import numpy.random as random
import scanpy as sc

large_root = "/home/breanne/sourcedata/limb_bud"

# read the data from 10x .mtx:
print("Reading Data")

new_data = sc.read_10x_mtx(large_root + "/e13")
newframe= pd.DataFrame(data=new_data.X.todense(), index=new_data.obs_names, columns=new_data.var_names)
newframe = newframe.T
newframe.to_csv(large_root + "/e13_limbbud.csv")

new_data2 = sc.read_10x_mtx(large_root + "/e12")
newframe2= pd.DataFrame(data=new_data2.X.todense(), index=new_data2.obs_names, columns=new_data2.var_names)
newframe2 = newframe2.T
newframe2.to_csv(large_root + "/e12_limbbud.csv")

e11_data = sc.read_10x_mtx(large_root + "/e11")
new_frame3= pd.DataFrame(data=e11_data.X.todense(), index=e11_data.obs_names, columns=e11_data.var_names)
new_frame3 = new_frame3.T
new_frame3.to_csv(large_root + "/e11_limbbud.csv")

new_data4 = sc.read_10x_mtx(large_root + "/e13_5")
new_frame4= pd.DataFrame(data=new_data4.X.todense(), index=new_data4.obs_names, columns=new_data4.var_names)
new_frame4 = new_frame4.T
new_frame4.to_csv(large_root + "/e13_5_limbbud.csv")

clust_list = [newframe, newframe2, new_frame3, new_frame4]
label = ['e13_limbbud','e12_limbbud','e11_limbbud','e13_5_limbbud']

for d, l in zip(clust_list, label):
    gene_headers = np.array(d.index)
    print("Cluster:" + l)
    data = d.to_numpy(dtype=float)
    avg_oall = np.mean(data, 1)
    coeffv0s = variation(data, axis=1)
    data[data == 0] = np.nan

    #avg_oall[avg_oall == 0] = np.nan
    avg_level = np.nanmean(data, 1)
    count_depth = np.sum(data,1) #how many times each mRNA is counted across all cells
    num_cells = np.count_nonzero(~np.isnan(data),1) #how many cells each mRNA is detected in
    NT = data.shape[1] #total ncells there are 17k genes in hydra! 25k cells

    Pc = 0.05
    maxv = (np.nanmax(avg_oall))
    minv = (np.nanmin(avg_oall))
    ma = np.logspace(np.log10(minv), np.log10(maxv), 5000)
    ms = [x / Pc for x in ma]
    Yexp = NT*(1 - (np.power((1 - Pc), ms)))

    n_starting_avg = np.round(avg_oall*NT/Pc)
    n_i_starting = np.round(n_starting_avg/NT)
    Pr_zero = (1-Pc)**n_i_starting

    numc = str(NT)
    getp = np.vstack((gene_headers, num_cells, avg_oall))
    getpar = np.transpose(getp)
    getparams = pd.DataFrame(getpar)
    getparams.to_csv(large_root + "/" + l + "gene_params.csv", header=None, index=None)
    base_fn = (l + numc + "gene_params.txt")
    with open(os.path.join(large_root, base_fn),'w') as outfile:
        getparams.to_string(outfile, header=None, index=None)

