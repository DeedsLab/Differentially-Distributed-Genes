import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.colors
import scipy.stats as ss
import seaborn as sns
from kneed import KneeLocator
from scipy.stats import variation
import scipy.stats as st


large_root = r"/media/timothyhamilton/My Passport1/Tim_Hamilton/Hydra/subclusters"
print("Reading data...")
label='LV5_clusters_hydra'
#all data
#cytT = pd.read_csv(large_root + r"/Zheng_FACS_pbmcs/cytotoxic_t.csv", index_col=0, usecols = [0,1])
#zheng = pd.read_csv(large_root + r"/Zheng_FACS_pbmcs/zheng9_5k/zheng5k_dropzeros.csv", index_col=0,usecols = [0,1])
#hydra = pd.read_csv(large_root + r"/R01_Ref_Datasets/Hydra_Results/Hydra_Raw_Cell_Types.csv", nrows=1)
#hydra = hydra.T
#bladder = pd.read_csv(large_root + r"/Mouse_bladder/mouse_bladder.csv", index_col=0,usecols = [0,1]) #gene by cell
#kidney = pd.read_csv(large_root + r"/Mouse_kidney/mkidney.csv", index_col=0,usecols = [0,1])
#planaria = pd.read_csv(large_root + r"/Planarian_atlas/planarian.csv", index_col=0,usecols = [0,1]) #nrows=1
#ddgs

c0 = pd.read_csv(large_root + r"/subcluster_0.csv", index_col=0)
c1 = pd.read_csv(large_root + r"/subcluster_1.csv", index_col=0)
c2 = pd.read_csv(large_root + r"/subcluster_2.csv", index_col=0)
c3 = pd.read_csv(large_root + r"/subcluster_3.csv", index_col=0)
c4 = pd.read_csv(large_root + r"/subcluster_4.csv", index_col=0)
c5 = pd.read_csv(large_root + r"/subcluster_5.csv", index_col=0)
c6 = pd.read_csv(large_root + r"/subcluster_6.csv", index_col=0)
c7 = pd.read_csv(large_root + r"/subcluster_7.csv", index_col=0)
c8 = pd.read_csv(large_root + r"/subcluster_8.csv", index_col=0)

all_list = [c0,c1,c2,c3,c4,c5,c6,c7,c8]

logcount = []
logfrac = []
for data in all_list:
    dropzero = data[np.all(data == 0, axis=1)].index
    data = data.drop(dropzero, 'index')  # drop genes that are not detected in any cells
    totalcells = data.shape[1]
    avgc = np.log10(np.mean(data,1))
    data[data == 0] = np.nan
    num_cells = data.count(axis=1)  # number of cells each mRNA is detected in
    num_cells[np.isnan(num_cells)] = 0
    frac_cells = np.log10(num_cells / totalcells)
    avgc[np.isnan(avgc)] = 0
    logcount.append(avgc)
    logfrac.append(frac_cells)

fs = 12
fig = plt.figure(figsize=(5, 8))
ax = fig.add_subplot(111)
#ax = sns.boxplot(data=[logcount[0],logcount[1],logcount[2],logcount[3],logcount[4],logcount[5],logcount[6],logcount[7],logcount[8]], width = 0.4, fliersize = 2)
ax = sns.violinplot(data=[logcount[0],logcount[1],logcount[2],logcount[3],logcount[4],logcount[5],logcount[6],logcount[7],logcount[8]], width = 1)
#plt.setp(ax1.collections)
ax.set_title("ncounts per cluster", fontsize=fs)
ax.set_ylabel("freq")
ax.set_xlabel("cluster")
#ax.set_yscale('log')
#ax.set_ylim(0, 0.6)
#plt.grid()
plt.savefig(large_root + "/" + label + "fig1_ncounts.png")
plt.savefig(large_root + "/" + label + "fig1_ncounts.eps")

fig = plt.figure(figsize=(5, 8))
ax = fig.add_subplot(111)
#ax = sns.boxplot(data=[logfrac[0],logfrac[1],logfrac[2],logfrac[3],logfrac[4],logfrac[5],logfrac[6],logfrac[7],logfrac[8]], width = 0.4, fliersize = 2)
ax = sns.violinplot(data=[logfrac[0],logfrac[1],logfrac[2],logfrac[3],logfrac[4],logfrac[5],logfrac[6],logfrac[7],logfrac[8]], width = 1)
#plt.setp(ax1.collections)
ax.set_title("frac cells per cluster", fontsize=fs)
ax.set_ylabel("freq")
ax.set_xlabel("cluster")
#ax.set_yscale('log')
#ax.set_ylim(0, 0.6)
#plt.grid()
plt.savefig(large_root + "/" + label + "fig1_fraccells.png")
plt.savefig(large_root + "/" + label + "fig1_fraccells.eps")