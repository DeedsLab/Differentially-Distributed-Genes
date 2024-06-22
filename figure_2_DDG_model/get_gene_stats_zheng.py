import pandas as pd
import numpy as np
import os as os
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy.stats import variation
from scipy.stats import binom
import numpy.random as random

large_root = r"C:\Users\bsparta\Desktop\data\clustering_validation\raw_data\zheng"

bcells = pd.read_csv(large_root + r"\b_cells.csv")
bcells = bcells.set_index('Unnamed: 0')
bcells = bcells.T #must be transposed

thelp = pd.read_csv(large_root + r"\cd4_t_helper.csv")
thelp = thelp.set_index('Unnamed: 0')
thelp = thelp.T #must be transposed

cd34 = pd.read_csv(large_root + r"\cd34.csv")
cd34 = cd34.set_index('Unnamed: 0')
cd34 = cd34.T #must be transposed

cd56_nk = pd.read_csv(large_root + r"\cd56_nk.csv")
cd56_nk = cd56_nk.set_index('Unnamed: 0')
cd56_nk = cd56_nk.T #must be transposed

cytotoxic_t = pd.read_csv(large_root + r"\cytotoxic_t.csv")
cytotoxic_t = cytotoxic_t.set_index('Unnamed: 0')
cytotoxic_t = cytotoxic_t.T #must be

ercc_filtered = pd.read_csv(large_root + r"\ercc_filtered.csv")
ercc_filtered = ercc_filtered.set_index('Unnamed: 0')
ercc_filtered = ercc_filtered.T #must be

memory_t = pd.read_csv(large_root + r"\memory_t.csv")
memory_t = memory_t.set_index('Unnamed: 0')
memory_t = memory_t.T #must be

monocytes = pd.read_csv(large_root + r"\monocytes.csv")
monocytes = monocytes.set_index('Unnamed: 0')
monocytes = monocytes.T #must be

naive_cytotoxic = pd.read_csv(large_root + r"\naive_cytotoxic.csv")
naive_cytotoxic = naive_cytotoxic.set_index('Unnamed: 0')
naive_cytotoxic = naive_cytotoxic.T #must be

naive_t = pd.read_csv(large_root + r"\naive_t.csv")
naive_t = naive_t.set_index('Unnamed: 0')
naive_t = naive_t.T #must be

regt = pd.read_csv(large_root + r"\regulatory_t.csv")
regt = regt.set_index('Unnamed: 0')
regt = regt.T #must be transposed

clust_list = [bcells, thelp, cd34, cd56_nk, cytotoxic_t, ercc_filtered,
              memory_t, monocytes, naive_cytotoxic, naive_t, regt]
label = ['b_cells','cd4_t_helper', 'cd34','cd56_nk', 'cytotoxic_t','ercc_filtered',
         'memory_t','monocytes','naive_cytotoxic', 'naive_t','regulatory_t']



for d, l in zip(clust_list, label):
    gene_headers = np.array(d.index)
    print("Cluster:" + l)
    data = d.to_numpy(dtype=float)
    avg_oall = np.mean(data, 1)
    coeffv0s = variation(data, axis=1)
    data[data == 0] = np.nan
    avg_oall[avg_oall == 0] = np.nan
    avg_level = np.nanmean(data, 1)
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
    #getparams.to_csv(large_root + "/" + l + "gene_params.csv", header=None, index=None)
    base_fn = (l + numc + "gene_params.txt")
    with open(os.path.join(large_root, base_fn),'w') as outfile:
        getparams.to_string(outfile, header=None, index=None)

    count_max = np.nanmax(data,1) #max across cells
    coeffv = variation(data, axis=1, nan_policy='omit')

    cm = plt.cm.get_cmap('cubehelix')
    cm2 = plt.cm.get_cmap('viridis')
    cm3 = plt.cm.get_cmap('cividis')
    cm4 = plt.cm.get_cmap('plasma')

    #heatmap
    fs = 30
    fig = plt.figure(figsize = (20,20))
    ax = fig.add_subplot(111)
    ax.set_title("mRNA detection level in" + l, fontsize=fs)
    ax.set_ylabel("genes")
    ax.set_xlabel("cells")
    hm = ax.pcolormesh(data, cmap=cm2, vmin=0, vmax = 1500)
    plt.colorbar(hm, extend='both')
    plt.savefig(large_root + "/" + l + "gc_heatmap.png")

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
