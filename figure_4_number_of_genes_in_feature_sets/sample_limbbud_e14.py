import pandas as pd
import numpy as np
import os as os
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy.stats import variation
from scipy.stats import binom
import numpy.random as random
import math
import scipy.stats as ss

large_root = r"/home/breanne/sourcedata/limb_bud"
#large_root = r"C:\Users\Breanne Sparta\Desktop\from posieden"
#test = pd.read_csv(large_root + r"/e12_limbbuddrop0.csv", index_col=0, usecols = [1,2,3])
rep_list = [1,2,3,4,5,6,7,8,9,10]
#rep_list = [1]

#e13 = pd.read_csv(large_root + r"/e13_limbbuddrop0.csv", index_col=0) #gene by cell
#e12 = pd.read_csv(large_root + r"/e12_limbbuddrop0.csv", index_col=0)
#e11 = pd.read_csv(large_root + r"/e11_limbbuddrop0.csv", index_col=0)
#e13_5 = pd.read_csv(large_root + r"/e13_5_limbbuddrop0.csv", index_col=0)
e14 = pd.read_csv(large_root + r"/e14_limbbuddrop0.csv", index_col=0)
#e13b = pd.read_csv(large_root + r"/e13b_limbbuddrop0.csv", index_col=0)
#e15 = pd.read_csv(large_root + r"/e15_limbbuddrop0.csv", index_col=0)
#e10_5 = pd.read_csv(large_root + r"/e10_5_limbbuddrop0.csv", index_col=0)

#datalist = [e13,e12,e11,e13_5,e14,e13b,e15,e10_5]
#labels = ["e13","e12","e11","e13_5","e14","e13b","e15","e10_5"]

data = e14
label = '_e14_'
cellheads = data.columns.tolist()
#first filter bulk of outliers that would affect median filter
genes_pcell = np.count_nonzero(data,0)
counts_pcell = np.sum(data, 0)  # reads per cell
filt = pd.DataFrame(index = cellheads)
filt['genes'] = genes_pcell
filt['counts'] = counts_pcell
# could filter outliers by percentile within window of 5 genes?
# could aso do a ratio
# or fit a polynomial and just add a value above that that is cut off
droplist = filt.index[(filt['genes']<=6000) & (filt['counts']>45000)].tolist()
droplist1 = filt.index[(filt['genes']<=5100) & (filt['counts']>39000)].tolist()
droplist2 = filt.index[(filt['genes']<=3900) & (filt['counts']>27000)].tolist()
droplist3 = filt.index[(filt['genes']<=2900) & (filt['counts']>19500)].tolist()
droplist4 = filt.index[(filt['genes']<=2100) & (filt['counts']>12000)].tolist()
droplist5 = filt.index[(filt['genes']<=1300) & (filt['counts']>9000)].tolist()
droplist6 = filt.index[(filt['genes']<=900) & (filt['counts']>7500)].tolist()
droplist7 = filt.index[(filt['genes']<=500) & (filt['counts']>1750)].tolist()

for item in droplist1:
    droplist.append(item)
for item in droplist2:
    droplist.append(item)
for item in droplist3:
    droplist.append(item)
for item in droplist4:
    droplist.append(item)
for item in droplist5:
    droplist.append(item)
for item in droplist6:
    droplist.append(item)
for item in droplist7:
    droplist.append(item)

data = data.drop(droplist, 'columns')
filt = filt.drop(droplist, 'index')

data.to_csv(large_root + "/resample/" + label + "_filt.csv")

datalist = [data]
labels = ['_e14_']

fs = 14
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(121)
ax.set_title(label, fontsize=fs)
ax.set_ylabel("counts per cell before")
ax.set_xlabel("genes per cell before")
sc = plt.scatter(genes_pcell, counts_pcell, s=30, marker="o")
plt.tight_layout()
plt.grid(markevery = 100)

ax = fig.add_subplot(122)
ax.set_title(label, fontsize=fs)
ax.set_ylabel("counts per cell after")
ax.set_xlabel("genes per cell after")
sc = plt.scatter(filt.genes, filt.counts, s=30, marker="o")
plt.tight_layout()
plt.grid()

plt.savefig(large_root + "/resample/" + label + "countfilter.png")
plt.savefig(large_root + "/resample/" + label + "countfilter.eps")
fig.clear()

#need to apply 2d kernal filtering to remove weird outlier cells
#or medfilt2d scipy median filter 2d array
    #bins of 1000, remove 8 outliers?

for group,l in zip(datalist, labels):
    for rep in rep_list:
        rep = str(rep)
        s=len(group.columns)
        randnums= np.random.randint(1,s,5000)
        d = group.iloc[:,randnums]
        dropzero = d[np.all(d == 0, axis=1)].index
        d2 = d.drop(dropzero, 'index')  # drop genes that are not detected in any cells
        #d = d.T #must be transposed
        #d2.to_csv(large_root + "/resample" + "_rep_" + rep + l + "limb_bud.csv")

        gene_headers = np.array(d2.index)
        print("Cluster:" + l)
        data = d2.to_numpy(dtype=float)
        avg_oall = np.mean(data, 1)
        coeffv0s = variation(data, axis=1)
        data[data == 0] = np.nan

        #avg_oall[avg_oall == 0] = np.nan
        avg_level = np.nanmean(data, 1)
        count_depth = np.sum(data,1) #how many times each mRNA is counted across all cells
        num_cells = np.count_nonzero(~np.isnan(data),1) #how many cells each mRNA is detected in
        NT = s #total ncells there are 17k genes in hydra! 25k cells

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
        getparams.to_csv(large_root + "/resample/" + l + "gene_params.csv", header=None, index=None)
        base_fn = ("rep" + rep + l + numc + "gene_params.txt")
        with open(os.path.join(large_root, base_fn),'w') as outfile:
            getparams.to_string(outfile, header=None, index=None)

        count_max = np.nanmax(data,1) #max across cells
        coeffv = variation(data, axis=1, nan_policy='omit')

        #### log plts

        cm = plt.cm.get_cmap('cubehelix')
        fs = 10
        fig = plt.figure(figsize=(5, 5))
        ax = fig.add_subplot(111)
        ax.set_title(l, fontsize=fs)
        ax.set_ylabel("number of cells")
        ax.set_xlabel("average mRNA over all cells")
        ax.loglog()
        z = coeffv0s
        sc2 = plt.scatter(avg_oall, num_cells, c=z, vmin=0, vmax=np.nanmax(z) + (.15 * np.nanmax(z)), s=30, edgecolors='none', marker=".", cmap=cm)
        plt.colorbar(sc2, label="coeff variation", orientation='horizontal', extend='both')
        sc2 = plt.plot(ma, Yexp, linewidth=3)

        plt.savefig(large_root + "/resample/_rep_" + rep + l + "Log countd numg scatter.png")
        plt.savefig(large_root + "/resample/_rep_" + rep + l + "Log countd numg scatter.eps")
        fig.clear()
