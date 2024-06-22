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

def rename(A):
    # for gene list A return gene names with spaces instead of - or . characters
    newlist = []
    for gene in A:
        newgene = str(gene)
        #a = newgene.replace('.', '').replace('-', '') #zheng has a gene where if you remove both . and -, two genes have the same name
        # probably where nan comes from in hvg data
        a = newgene.replace('.', '-') #R changes - to .
        newlist.append(a)
    all_headers = set(np.array(newlist))
    return (all_headers)

#calc jaccard index
def WC_similarity(list1):
    s1 = set(list1)
    s2 = set(gene_list[1]) #this is WC list
    return float(len(s1.intersection(s2)) / len(s1.union(s2)))

large_root = "/media/timothyhamilton/data1/Tim_Hamilton/Zheng_FACS_pbmcs/zheng9_5k"
label = 'zheng5k'

total_data = pd.read_csv("/media/timothyhamilton/data1/Tim_Hamilton/Zheng_FACS_pbmcs/zheng9_5k/zheng5k_dropzeros.csv", index_col=0,) #gene x cells
WC = pd.read_csv("/media/timothyhamilton/data1/Tim_Hamilton/Zheng_FACS_pbmcs/zheng9_5k/zheng9_5k_WC.csv", index_col=0, usecols = [0,1]) #g x c
HVG = pd.read_csv("/media/timothyhamilton/data1/Tim_Hamilton/Zheng_FACS_pbmcs/zheng9_5k/HVGs.csv",nrows=1)
DDG = pd.read_csv("/media/timothyhamilton/data1/Tim_Hamilton/Zheng_FACS_pbmcs/zheng9_5k/zheng9_5k_DDGs.csv", index_col=0, usecols = [0,1])# g x c
townes = pd.read_csv("/media/timothyhamilton/data1/Tim_Hamilton/Zheng_FACS_pbmcs/zheng9_5k/zheng5k_dropzeros_townes.csv",index_col = 0,nrows=1)
mdrop = pd.read_csv("/media/timothyhamilton/data1/Tim_Hamilton/Zheng_FACS_pbmcs/zheng9_5k/zheng5k_dropzeros_M3Drop.csv",index_col = 0,nrows=1)
WC = WC.T
DDG = DDG.T
total_data = total_data.T

# make gene sets
all_set= set(np.array(total_data.columns))
DDG_set= set(np.array(DDG.columns))
WC_set= set(np.array(WC.columns))
M3Drop_set= set(np.array(mdrop.columns))
Townes_set= set(np.array(townes.columns))
HVG_set = set(np.array(HVG.columns))

#format data
title_list = ["all_genes","WC","DDGs","HVGs","M3Drop","Townes"]
val_list = [all_set, WC_set, DDG_set, HVG_set, M3Drop_set, Townes_set]
gene_list = [rename(A) for A in val_list]

#j = len(gene_list[1])
#print("formatted gene list wc"+str(j))

cell_headers = np.array(total_data.index)
total_data = total_data.values
total_data= pd.DataFrame(total_data, index = cell_headers, columns = gene_list[0])
title_list2 = title_list[1:]
gene_list2 = gene_list[1:]

#calc jaccard index with WC group
jindex = [WC_similarity(list1) for list1 in gene_list2]
fig = plt.figure(figsize=(12, 8))
ax1 = fig.add_subplot(111)
jindex_frame = pd.DataFrame({'gene_set': title_list2, 'AJI': jindex})
jindex_frame.to_csv(large_root + "/" + label + "jaccardindex_w_WCs.csv")
fs = 12
ax1 = jindex_frame.plot.bar(x='gene_set', y='AJI', rot=0, width = 0.75)
ax1.set_ylabel("Jaccard Index")
ax1.set_xlabel("set compared to WC genes")
ax1.set_title("proportional overlap with WC genes")
ax1.set_ylim(0,1)
plt.grid()
plt.savefig(large_root + "/" + label + "jaccardI_w_WC.eps")
plt.savefig(large_root + "/" + label + "jaccardI_w_WC.png")

#venn across all groups
val_dict= dict(list(zip(title_list2,gene_list2)))
fig = plt.figure(figsize=(15, 8))
ax1 = fig.add_subplot(111)
ax1.set_title("Overlap of Gene_Sets")
venn(val_dict, cmap="viridis", fontsize=12, legend_loc="upper left", ax=ax1)
ax1.set_rasterized(True)
plt.savefig(large_root + "/" + label + "gene_overlap.eps")
plt.savefig(large_root + "/" + label + "gene_overlap.png")

print("files updated")
data_list = [total_data]
cellname = ['_all_zheng5k']

for cell_group, clabel in zip(data_list,cellname):
    #WC_data = cell_group.filter(WC_headers, axis=0)
    data2 = cell_group.T #gene x cell
    avg_oall = np.mean(data2, 1) #avg mRNA exp across all cells
    data2[data2 == 0] = np.nan
    num_cells = data2.count(axis=1)  # number of cells each mRNA is detected in
    NT = data2.shape[1]  # total ncells there are 17k genes in hydra! 25k cells
    n_cells = str(NT)
    data2 = data2.fillna(0)

    Pc = 0.05
    avg_oall[avg_oall == 0] = np.nan
    maxv = (np.nanmax(avg_oall))
    minv = (np.nanmin(avg_oall))
    ma = np.logspace(np.log10(minv), np.log10(maxv), 5000)
    ms = [x / Pc for x in ma]
    Yexp = NT * (1 - (np.power((1 - Pc), ms)))
    n_starting_avg = np.round(avg_oall * NT / Pc)
    n_i_starting = np.round(n_starting_avg / NT)
    Pr_zero = (1 - Pc) ** n_i_starting

    # here, plot all genes and then highlight gene list genes
    #then also do an overlap one
    for gene_group, glabel in zip(gene_list2, title_list2):
        # DDG set overlaps
        fig = plt.figure(figsize=(15, 8))
        ax = fig.add_subplot(111)
        venn2([set(gene_list[2]), set(gene_group)], set_labels=('DDGs', glabel))
        ax.set_rasterized(True)
        plt.savefig(large_root + "/" + label + "gene_overlap_DDG" + glabel + ".png")
        plt.savefig(large_root + "/" + label + "gene_overlap_DDG" + glabel + ".eps")

        # WC set overlaps
        fig = plt.figure(figsize=(15, 8))
        ax = fig.add_subplot(111)
        venn2([set(gene_list[1]), set(gene_group)], set_labels=('WC', glabel))
        ax.set_rasterized(True)
        plt.savefig(large_root + "/" + label + "gene_overlap_WC" + glabel + ".png")
        plt.savefig(large_root + "/" + label + "gene_overlap_WC" + glabel + ".eps")

        #scatter gene sets with intersections
        test_data = data2.filter(gene_group, axis=0)
        avg_oall2 = np.mean(test_data, 1)  # avg mRNA exp across all cells
        test_data[test_data == 0] = np.nan
        num_cells2 = test_data.count(axis=1)  # number of cells each mRNA is detected in

        #intersect w WC
        test_data2 = data2.filter(gene_list[1], axis=0) #hardcode which gene set to compare
        avg_oall3 = np.mean(test_data2, 1)  # avg mRNA exp across all cells
        test_data2[test_data2 == 0] = np.nan
        num_cells3 = test_data2.count(axis=1)  # number of cells each mRNA is detected in

        intersectWC = gene_list[1].intersection(gene_group)
        i = len(intersectWC)
        print("intersect"+str(i))
        intersectWC = list(intersectWC)
        inter_data = data2.filter(intersectWC, axis=0)
        avg_oall4 = np.mean(inter_data, 1)  # avg mRNA exp across all cells
        inter_data[inter_data == 0] = np.nan
        num_cells4 = inter_data.count(axis=1)  # number of cells each mRNA is detected in

        fs = 14
        fig = plt.figure(figsize=(12, 6))
        ax = fig.add_subplot(121)
        ax.set_title(label + glabel, fontsize=fs)
        ax.set_ylabel("number of cells")
        ax.set_xlabel("average mRNA count over all cells")
        ax.loglog()
        sc1 = plt.scatter(avg_oall, num_cells, c='black', s=80, edgecolors='None',
                          marker=".", label='all genes')
        sc1 = plt.scatter(avg_oall2, num_cells2, c='m', s=80, edgecolors='None',
                          marker=".", label=glabel)
        ax.legend()
        sc1 = plt.plot(ma, Yexp, linewidth=3)
        plt.tight_layout()

        ax = fig.add_subplot(122)
        ax.set_title(label + glabel + "overlap with WC set", fontsize=fs)
        ax.set_ylabel("number of cells")
        ax.set_xlabel("average mRNA count over all cells")
        ax.loglog()
        sc1 = plt.scatter(avg_oall, num_cells, c='black', s=80, edgecolors='None',
                          marker=".", label = 'all genes')
        sc1 = plt.scatter(avg_oall2, num_cells2, c='darkgoldenrod', s=80, edgecolors='None',
                          marker=".", label = glabel)
        sc1 = plt.scatter(avg_oall3, num_cells3, c='blue', s=80, edgecolors='None',
                          marker=".", label= 'WC')
        sc1 = plt.scatter(avg_oall4, num_cells4, c='chartreuse', s=80, edgecolors='None',
                          marker=".", label='both')
        ax.legend()
        sc1 = plt.plot(ma, Yexp, linewidth=3)
        plt.tight_layout()
        plt.savefig(large_root + "/" + label + glabel + "scatter_set_intersect_WC.png")
        plt.savefig(large_root + "/" + label + glabel + "scatter_set_intersect_WC.eps")
        fig.clear()

    #make split violin plot
    #for each gene group plot avg mRNA count vs percent of cells its detected in
    #accumulate data into pandas frame
    #cat mRNA count and percent of cells into one column, add second column as variable to show which is which

avg_list = []
perc_list = []
c = 0

fig = plt.figure(figsize=(18, 6))
for gene_group, glabel in zip(gene_list, title_list):
    c = c+1
    #make list of dataframes
    test_data = data2.filter(gene_group, axis=0)
    avg_oall = np.mean(test_data, 1)  # avg mRNA exp across all cells
    test_data[test_data == 0] = np.nan
    num_cells = test_data.count(axis=1)  # number of cells each mRNA is detected in
    perc_cells = [(x / NT) * 100 for x in num_cells]

    count_frame = pd.DataFrame({'avg_count': avg_oall, 'perc_cells': perc_cells})
    count_data = count_frame.melt().assign(x=glabel)

    ax = fig.add_subplot(1,len(title_list),c)
    ax.set_title(label + glabel , fontsize=fs)
    ax.set_ylabel("geneset")
    ax.set_xlabel("value")
    #plt.ylim((10 ** -5, 10 ** 2))
    ax.set(yscale="log")
    ax = sns.violinplot(data=count_data, x='x', y='value',
                   hue='variable', split=True, inner='quart', scale = "count")
    plt.tight_layout()

    glabel = avg_oall
    avg_list.append(glabel)
    glabel = perc_cells
    perc_list.append(glabel)

plt.savefig(large_root + "/" + label + "count_stats.png")
plt.savefig(large_root + "/" + label + "count_stats.eps")
fig.clear()

fig = plt.figure(figsize=(18, 6))
ax = fig.add_subplot(121)
ax.set_title("average mRNA across cells for:" + label + "set", fontsize=fs)
ax.set_ylabel("average mRNA count")
ax.set_xlabel("geneset")
#ax.set(yscale="log")
ax = sns.violinplot(data=avg_list, scale = "count")
ax.set(xticklabels=title_list)

ax = fig.add_subplot(122)
ax.set_title("percent cells mRNA is detected in for:" + label + "set", fontsize=fs)
ax.set_ylabel("percent cell count per mRNA")
ax.set_xlabel("geneset")
#ax.set(yscale="log")
ax = sns.violinplot(data=perc_list, scale = "count")
ax.set(xticklabels=title_list)
plt.savefig(large_root + "/" + label + "compare_count_statsnl.png")
plt.savefig(large_root + "/" + label + "compare_count_statsnl.eps")