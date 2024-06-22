import pandas as pd
import numpy as np
import os as os
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy.stats import variation
from scipy.stats import binom
import numpy.random as random

large_root = r"C:\Users\bsparta\Desktop\data\clustering_validation\raw_data\zheng"

randnums= np.random.randint(1,10001,5000)
bcells = pd.read_csv(large_root + r"\b_cells.csv", index_col = 0)
bcells = bcells.iloc[randnums,:]
bcells = bcells.T #must be transposed

randnums= np.random.randint(1,10001,5000)
thelp = pd.read_csv(large_root + r"\cd4_t_helper.csv", index_col = 0)
thelp = thelp.iloc[randnums,:]
thelp = thelp.T #must be transposed

randnums= np.random.randint(1,9232,5000)
cd34 = pd.read_csv(large_root + r"\cd34.csv", index_col = 0)
cd34 = cd34.iloc[randnums,:]
cd34 = cd34.T #must be transposed

randnums= np.random.randint(1,8385,5000)
cd56_nk = pd.read_csv(large_root + r"\cd56_nk.csv", index_col = 0)
cd56_nk = cd56_nk.iloc[randnums,:]
cd56_nk = cd56_nk.T #must be transposed

randnums= np.random.randint(1,10001,5000)
cytotoxic_t = pd.read_csv(large_root + r"\cytotoxic_t.csv", index_col = 0)
cytotoxic_t = cytotoxic_t.iloc[randnums,:]
cytotoxic_t = cytotoxic_t.T #must be

randnums= np.random.randint(1,10001,5000)
memory_t = pd.read_csv(large_root + r"\memory_t.csv", index_col = 0)
memory_t = memory_t.iloc[randnums,:]
memory_t = memory_t.T #must be

monocytes = pd.read_csv(large_root + r"\monocytes.csv", index_col = 0)
monocytes = monocytes.T #must be

randnums= np.random.randint(1,10001,5000)
naive_cytotoxic = pd.read_csv(large_root + r"\naive_cytotoxic.csv", index_col = 0)
naive_cytotoxic = naive_cytotoxic.iloc[randnums,:]
naive_cytotoxic = naive_cytotoxic.T #must be

randnums= np.random.randint(1,10001,5000)
naive_t = pd.read_csv(large_root + r"\naive_t.csv", index_col = 0)
naive_t = naive_t.iloc[randnums,:]
naive_t = naive_t.T #must be

randnums= np.random.randint(1,10001,5000)
regt = pd.read_csv(large_root + r"\regulatory_t.csv", index_col = 0)
regt = regt.iloc[randnums,:]
regt = regt.T #must be transposed

bcells.columns = ['b_cells' + str(col) for col in bcells.columns]
thelp.columns = ['cd4_t_helper' + str(col) for col in thelp.columns]
cd56_nk.columns = ['cd56_nk' + str(col) for col in cd56_nk.columns]
cytotoxic_t.columns = ['cytotoxic_t' + str(col) for col in cytotoxic_t.columns]
memory_t.columns = ['memory_t' + str(col) for col in memory_t.columns]
monocytes.columns = ['monocytes' + str(col) for col in monocytes.columns]
naive_cytotoxic.columns = ['naive_cytotoxic' + str(col) for col in naive_cytotoxic.columns]
naive_t.columns = ['naive_t' + str(col) for col in naive_t.columns]
regt.columns = ['regulatory_t' + str(col) for col in regt.columns]

bthelp = pd.merge(bcells, thelp, left_index=True, right_index=True,how='outer')
bthelp = bthelp.replace(np.nan,0)
#save
#consider making smaller mergedatafile, with 2k cells per group, do it
mergedata = pd.merge(bthelp, cd56_nk, left_index=True, right_index=True,how='outer')
mergedata = pd.merge(mergedata, regt, left_index=True, right_index=True,how='outer')
mergedata = pd.merge(mergedata, cytotoxic_t, left_index=True, right_index=True,how='outer')
mergedata = pd.merge(mergedata, memory_t, left_index=True, right_index=True,how='outer')
mergedata = pd.merge(mergedata, monocytes, left_index=True, right_index=True,how='outer')
mergedata = pd.merge(mergedata, naive_cytotoxic, left_index=True, right_index=True,how='outer')
mergedata = pd.merge(mergedata, naive_t, left_index=True, right_index=True,how='outer')
mergedata = mergedata.fillna(0)
mergedata.to_csv(large_root + "/zheng9_5k.csv")

clust_list = [mergedata]
label = ['zheng9_5k']

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

    count_max = np.nanmax(data,1) #max across cells
    coeffv = variation(data, axis=1, nan_policy='omit')

    cm = plt.cm.get_cmap('cubehelix')
    cm2 = plt.cm.get_cmap('viridis')
    cm3 = plt.cm.get_cmap('cividis')
    cm4 = plt.cm.get_cmap('plasma')


#### log plts
    fs = 10
    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(121)
    ax.set_title(l, fontsize=fs)
    ax.set_ylabel("number of cells")
    ax.set_xlabel("average mRNA in expressing cells")
    ax.loglog()
    z = coeffv
    sc = plt.scatter(avg_level, num_cells, c=z, s=30, vmin=0, vmax=np.nanmax(z) + (.15 * np.nanmax(z)), edgecolors='', marker=".", cmap=cm)
    plt.colorbar(sc, label='coeff variation in exp cells', orientation='horizontal', extend='both')
    sc = plt.plot(ma, Yexp, linewidth=3)

    ax = fig.add_subplot(122)
    ax.set_title(l, fontsize=fs)
    ax.set_ylabel("number of cells")
    ax.set_xlabel("average mRNA over all cells")
    ax.loglog()
    z = coeffv0s
    sc2 = plt.scatter(avg_oall, num_cells, c=z, vmin=0, vmax=np.nanmax(z) + (.15 * np.nanmax(z)), s=30, edgecolors='', marker=".", cmap=cm)
    plt.colorbar(sc2, label="coeff variation", orientation='horizontal', extend='both')
    sc2 = plt.plot(ma, Yexp, linewidth=3)

    plt.savefig(large_root + "/" + l + "Log countd numg scatter.png")
    plt.savefig(large_root + "/" + l + "Log countd numg scatter.eps")
    fig.clear()
