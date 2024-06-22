import pandas as pd
import numpy as np
import os as os
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy.stats import variation
from scipy.stats import binom
import numpy.random as random
import math

large_root = r"/home/breanne/sourcedata/Zheng9/cell_lines"
rep_list = [1,2,3,4,5,6,7,8,9,10]
#rep_list = [1]

intlist = [100,500,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000]

#cell x gene
bcells1 = pd.read_csv(large_root + r"/b_cells.csv", index_col=0)
thelp1 = pd.read_csv(large_root + r"/cd4_t_helper.csv", index_col=0)
cd341 = pd.read_csv(large_root + r"/cd34.csv", index_col=0)
cd56_nk1 = pd.read_csv(large_root + r"/cd56_nk.csv", index_col=0)
cytotoxic_t1 = pd.read_csv(large_root + r"/cytotoxic_t.csv", index_col=0)
memory_t1 = pd.read_csv(large_root + r"/memory_t.csv", index_col=0)
monocytes1 = pd.read_csv(large_root + r"/monocytes.csv", index_col=0)
naive_cytotoxic1 = pd.read_csv(large_root + r"/naive_cytotoxic.csv", index_col=0)
naive_t1 = pd.read_csv(large_root + r"/naive_t.csv", index_col=0)
regt1 = pd.read_csv(large_root + r"/regulatory_t.csv", index_col=0)

bcells1 = bcells1.T
thelp1 = thelp1.T
cd341 = cd341.T
cd56_nk1 = cd56_nk1.T
cytotoxic_t1 = cytotoxic_t1.T
memory_t1 = memory_t1.T
monocytes1 = monocytes1.T
naive_cytotoxic1 = naive_cytotoxic1.T
naive_t1 = naive_t1.T
regt1 = regt1.T

for rep in rep_list:
    for int in intlist:
        rep = str(rep)

        randnums= np.random.randint(1,10001,int)
        bcells = bcells1.iloc[:,randnums]

        randnums= np.random.randint(1,10001,int)
        thelp = thelp1.iloc[:,randnums]

        mmint = [math.floor(int*0.9232)]
        randnums= np.random.randint(1,9232,mmint)
        cd34 = cd341.iloc[:,randnums]

        mmmint = [math.floor(int * 0.8385)]
        randnums= np.random.randint(1,8385,mmmint)
        cd56_nk = cd56_nk1.iloc[:,randnums]

        randnums= np.random.randint(1,10001,int)
        cytotoxic_t = cytotoxic_t1.iloc[:,randnums]

        randnums= np.random.randint(1,10001,int)
        memory_t = memory_t1.iloc[:,randnums]

        ##scale monocytes
        mint = [math.floor(int*0.26)]
        randnums= np.random.randint(1,2600,mint)
        monocytes = monocytes1.iloc[:,randnums]

        randnums= np.random.randint(1,10001,int)
        naive_cytotoxic = naive_cytotoxic1.iloc[:,randnums]

        randnums= np.random.randint(1,10001,int)
        naive_t = naive_t1.iloc[:,randnums]

        randnums= np.random.randint(1,10001,int)
        regt = regt1.iloc[:,randnums]

        bcells.columns = ['b_cells' + str(col) for col in bcells.columns]
        thelp.columns = ['cd4_t_helper' + str(col) for col in thelp.columns]
        cd56_nk.columns = ['cd56_nk' + str(col) for col in cd56_nk.columns]
        cytotoxic_t.columns = ['cytotoxic_t' + str(col) for col in cytotoxic_t.columns]
        memory_t.columns = ['memory_t' + str(col) for col in memory_t.columns]
        monocytes.columns = ['monocytes' + str(col) for col in monocytes.columns]
        naive_cytotoxic.columns = ['naive_cytotoxic' + str(col) for col in naive_cytotoxic.columns]
        naive_t.columns = ['naive_t' + str(col) for col in naive_t.columns]
        regt.columns = ['regulatory_t' + str(col) for col in regt.columns]

        bcells = bcells.T
        thelp = thelp.T
        cd34 = cd34.T
        cd56_nk = cd56_nk.T
        cytotoxic_t = cytotoxic_t.T
        memory_t = memory_t.T
        monocytes = monocytes.T
        naive_cytotoxic = naive_cytotoxic.T
        naive_t = naive_t.T
        regt = regt.T

        comb = pd.concat([bcells, thelp, cd34, cd56_nk, cytotoxic_t, memory_t, monocytes, naive_cytotoxic, naive_t, regt])
        comb = comb.fillna(0)
        comb = comb.loc[:, (comb != 0).any()]  # drop genes that are not detected in any cell
        int2 = str(int)
        comb.to_csv(r"/media/timothyhamilton/Seagate Desktop Drive/Breanne/zheng_resampled" + int2 + "+rep_" + rep + "zheng9_resamp.csv")

        label = [int2 + "+rep_" + rep + "zheng_resamp"]
        clust_list = [mergedata]

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
            getparams.to_csv(large_root + "/resample" + l + "gene_params.csv", header=None, index=None)
            base_fn = (l + numc + "gene_params.txt")
            with open(os.path.join(large_root, base_fn),'w') as outfile:
                getparams.to_string(outfile, header=None, index=None)

            #count_max = np.nanmax(data,1) #max across cells
            #coeffv = variation(data, axis=1, nan_policy='omit')

        #### log plts

            #cm = plt.cm.get_cmap('cubehelix')
            #fs = 10
            #fig = plt.figure(figsize=(5, 5))
            #ax = fig.add_subplot(111)
            #ax.set_title(l, fontsize=fs)
            #ax.set_ylabel("number of cells")
            #ax.set_xlabel("average mRNA over all cells")
            #ax.loglog()
            #z = coeffv0s
            #sc2 = plt.scatter(avg_oall, num_cells, c=z, vmin=0, vmax=np.nanmax(z) + (.15 * np.nanmax(z)), s=30, edgecolors='none', marker=".", cmap=cm)
            #plt.colorbar(sc2, label="coeff variation", orientation='horizontal', extend='both')
            #sc2 = plt.plot(ma, Yexp, linewidth=3)

            #plt.savefig(large_root + "/resample" + l + "Log countd numg scatter.png")
            #plt.savefig(large_root + "/resample" + l + "Log countd numg scatter.eps")
            #fig.clear()
#mergedata.to_csv(large_root + "/resample" + int2 + "+rep_" + rep + "zheng9_resamp.csv")
