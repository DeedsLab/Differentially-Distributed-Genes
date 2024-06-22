import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.colors
import scipy.stats as ss
import seaborn as sns
from kneed import KneeLocator
from scipy.stats import variation
from venn import venn
from matplotlib_venn import venn2, venn3


#get gene stats
def get_gene_stats(test_data):
    avg_oall = np.mean(test_data, 1)  # avg mRNA exp across all cells
    coeffs = variation(test_data, axis=1)
    vari = np.var(test_data, axis=1)
    avg00s = avg_oall
    avg00s[avg00s == 0] = np.nan
    varmean = [x / m for x, m in zip(vari, avg00s)]
    #varmean = [x / m for x, m in zip(vari, avg_oall)]
    rankmean = ss.rankdata(avg_oall, method='min')
    rankrat = ss.rankdata(varmean, method='min')
    test_data[test_data == 0] = np.nan
    num_cells = test_data.count(axis=1)  # number of cells each mRNA is detected in
    ranknum = ss.rankdata(num_cells, method='min')
    stat_frame = pd.DataFrame({'avg_oall': avg_oall, 'num_cells': num_cells, 'vari': vari, 'varmean': varmean, 'rankmean': rankmean, 'rankrat': rankrat, 'ranknum': ranknum})
    stat_frame = stat_frame.fillna(0)
    return(stat_frame)

#pipeline for gene set vis for data for which there is no WC set

large_root = r"/media/timothyhamilton/My Passport1/Tim_Hamilton/Zheng_FACS_pbmcs/zheng9_5k"
label = 'zheng5k'

c0 = pd.read_csv(large_root + "/subclusters/subcluster_0.csv", index_col=0) #gene x cells use cols
c1 = pd.read_csv(large_root + "/subclusters/subcluster_1.csv", index_col=0) #gene x cells use cols
c2 = pd.read_csv(large_root + "/subclusters/subcluster_2.csv", index_col=0) #gene x cells use cols
c3 = pd.read_csv(large_root + "/subclusters/subcluster_3.csv", index_col=0) #gene x cells use cols
title_list = ["c0","c1","c2","c3"]
gene_list = [c0, c1, c2, c3]
print("files updated")

def meancount(test_data):
    avg_oall = np.mean(test_data, 1)  # avg mRNA exp across all cells
    return avg_oall

NT = c0.shape[1]
def fraccell(test_data):
    num_cells = np.count_nonzero(test_data, axis=1)  # number of cells each mRNA is detected in
    perc_cells = [(x / NT) *10 for x in num_cells]
    return perc_cells

meangcount = [meancount(list1) for list1 in gene_list]
fraccells = [fraccell(list1) for list1 in gene_list]
gcount = [np.log(list1) for list1 in meangcount]
fcells = [np.log(list1) for list1 in fraccells]

fig = plt.figure(figsize=(18, 6))
ax = fig.add_subplot(121)
fs = 12
ax.set_title("average mRNA across cells for:" + label + "set", fontsize=fs)
ax.set_ylabel("average mRNA count")
ax.set_xlabel("gene cluster")
ax.set(yscale="log")
ax.set_ylim(0.00001,80)
ax = sns.violinplot(data=gcount)
ax.set(xticklabels=title_list)

ax = fig.add_subplot(122)
ax.set_title("percent cells mRNA is detected in for:" + label + "set", fontsize=fs)
ax.set_ylabel("percent cell count per mRNA")
ax.set_xlabel("gene cluster")
ax.set(yscale="log")
ax.set_ylim(0.001,100)
ax = sns.violinplot(data=fcells)
ax.set(xticklabels=title_list)
plt.savefig(large_root + "/" + label + "compare_count_stats_nl.png")
plt.savefig(large_root + "/" + label + "compare_count_stats_nl.eps")