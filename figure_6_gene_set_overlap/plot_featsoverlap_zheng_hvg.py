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
        a = newgene.replace('.', '-') #R changes - to .
        newlist.append(a)
    all_headers = np.array(newlist)
    return (all_headers)

#calc jaccard index
def WC_similarity(list1):
    s1 = set(list1)
    s2 = set(gene_list[1]) #this is DDG list in this case
    return float(len(s1.intersection(s2)) / len(s1.union(s2)))

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

large_root = r"/home/breannesparta/ddgs/zheng9"
label = 'zheng5k'

total_data = pd.read_csv("/home/breannesparta/ddgs/zheng9/cellgene/zheng5k_dropzeros.csv", index_col=0) #gene x cells use cols
DDG1 = pd.read_csv(large_root + "/cellgene/zheng9pcap1DDG.csv", index_col=0, nrows=1) #c x g
HVG2 = pd.read_csv(large_root + "/cellgene/zheng92k_hvgs.csv", index_col=0, nrows=1) #c x g
HVG8 = pd.read_csv(large_root + "/cellgene/zheng98k_hvgs.csv", index_col=0, nrows=1) #c x g
total_data = total_data.T

all_set= rename(total_data.columns.tolist())
print(all_set)
hvgs2k = rename(HVG2.columns.tolist())
hvgs8k = rename(HVG8.columns.tolist())
ddgsp1 = rename(DDG1.columns.tolist())

#HVG_set = [x for x in HVG_first if x in all_set]

#gene in second position will be reference list
title_list = ["all_genes","2k_HVGs","8k_HVGs","pcap1_DDGs"]
gene_list = [all_set, hvgs2k, hvgs8k, ddgsp1]

cell_headers = np.array(total_data.index)
total_data = total_data.values
total_data= pd.DataFrame(total_data, index = cell_headers, columns = gene_list[0])
title_list2 = title_list[1:]
gene_list2 = gene_list[1:]

"""#calc jaccard index with WC group
jindex = [WC_similarity(list1) for list1 in gene_list2]
fig = plt.figure(figsize=(12, 8))
ax1 = fig.add_subplot(111)
jindex_frame = pd.DataFrame({'gene_set': title_list2, 'AJI': jindex})
jindex_frame.to_csv(large_root + "/" + label + "jaccardindex_w_DDGs.csv")
fs = 12
ax1 = jindex_frame.plot.bar(x='gene_set', y='AJI', rot=0, width = 0.75)
ax1.set_ylabel("Jaccard Index")
ax1.set_xlabel("set compared to WC genes")
ax1.set_title("proportional overlap with WC genes")
ax1.set_ylim(0,1)
plt.grid()
plt.savefig(large_root + "/" + label + "jaccardI_w_WC.eps")
plt.savefig(large_root + "/" + label + "jaccardI_w_WC.png")"""


print("files updated")

data_list = [total_data]
cellname = ['_raw_counts_']

for cell_group, clabel in zip(data_list,cellname):
    data2 = cell_group.T #gene x cell
    test_data = data2
    stat_frame_all = get_gene_stats(test_data)

    avg_oall2 = stat_frame_all.avg_oall
    Pc = 0.05
    avg_oall2[avg_oall2 == 0] = np.nan
    maxv = (np.nanmax(avg_oall2))
    minv = (np.nanmin(avg_oall2))
    ma = np.logspace(np.log10(minv), np.log10(maxv), 5000)
    ms = [x / Pc for x in ma]
    NT = data2.shape[1]  # total ncells there are 17k genes in hydra! 25k cells
    Yexp = NT * (1 - (np.power((1 - Pc), ms)))
    n_starting_avg = np.round(avg_oall2 * NT / Pc)
    n_i_starting = np.round(n_starting_avg / NT)
    Pr_zero = (1 - Pc) ** n_i_starting

    # here, plot all genes and then highlight gene list genes
    #then also make scatter plot with intersect of sets
    for gene_group, glabel in zip(gene_list2, title_list2):
        # Ven diagram set overlaps
        fig = plt.figure(figsize=(15, 8))
        ax = fig.add_subplot(111)
        venn2([set(gene_list[1]), set(gene_group)], set_labels=('HVGs2k', glabel))
        ax.set_rasterized(True)
        plt.savefig(large_root + "/" + label + clabel+ "gene_overlap_hvg2k_and_" + glabel + ".png")
        plt.savefig(large_root + "/" + label + clabel+ "gene_overlap_hvg2k_and_" + glabel + ".eps")
        fig.clear()

        """"# WC set overlaps
        fig = plt.figure(figsize=(15, 8))
        ax = fig.add_subplot(111)
        venn2([set(gene_list[1]), set(gene_group)], set_labels=('WC', glabel))
        ax.set_rasterized(True)
        plt.savefig(large_root + "/" + label + "gene_overlap_WC" + glabel + ".png")
        plt.savefig(large_root + "/" + label + "gene_overlap_WC" + glabel + ".eps")
        fig.clear()"""

        #scatter gene sets with intersections
        test_group = pd.Series(list(gene_group))
        stat_frame_gg = stat_frame_all.loc[test_group, :]

        #intersect w test group in this case hvg 2k
        ddg_group = pd.Series(list(gene_list[1]))
        stat_frame_ddgs = stat_frame_all.loc[ddg_group,:]

        set_2 = frozenset(gene_list[1])
        intersectWC = [x for x in gene_group if x in set_2]
        i = len(intersectWC)
        print("intersect"+str(i))
        intersectWC = list(intersectWC)
        stat_frame_int = stat_frame_all.filter(intersectWC, axis=0)

        fs = 14
        fig = plt.figure(figsize=(12, 6))

        ax = fig.add_subplot(121)
        ax.set_title(label + glabel, fontsize=fs)
        ax.set_ylabel("number of cells")
        ax.set_xlabel("average mRNA count over all cells")
        ax.loglog()
        sc1 = plt.scatter(stat_frame_all.avg_oall, stat_frame_all.num_cells, c='black', s=80, edgecolors='None',
                          marker=".", label = 'all genes')
        sc1 = plt.scatter(stat_frame_gg.avg_oall, stat_frame_gg.num_cells, c='deeppink', s=80, edgecolors='None',
                          marker=".", label = glabel)
        ax.legend()
        sc1 = plt.plot(ma, Yexp, linewidth=3)
        plt.tight_layout()

        ax = fig.add_subplot(122)
        ax.set_title(label + clabel + glabel + "overlap with WC set", fontsize=fs)
        ax.set_ylabel("number of cells")
        ax.set_xlabel("average mRNA count over all cells")
        ax.loglog()
        sc1 = plt.scatter(stat_frame_all.avg_oall, stat_frame_all.num_cells, c='black', s=80, edgecolors='None',
                          marker=".", label = 'all genes')
        sc1 = plt.scatter(stat_frame_gg.avg_oall, stat_frame_gg.num_cells, c='deeppink', s=80, edgecolors='None',
                          marker=".", label = glabel)
        sc1 = plt.scatter(stat_frame_ddgs.avg_oall, stat_frame_ddgs.num_cells, c='deepskyblue', s=80, edgecolors='None',
                          marker=".", label= title_list2[1])
        sc1 = plt.scatter(stat_frame_int.avg_oall, stat_frame_int.num_cells, c='gold', s=80, edgecolors='None',
                          marker=".", label='both')
        ax.legend()
        sc1 = plt.plot(ma, Yexp, linewidth=3)
        plt.tight_layout()
        plt.savefig(large_root + "/" + label + glabel + clabel+ "scatter_set_intersect.png")
        plt.savefig(large_root + "/" + label + glabel + clabel+ "scatter_set_intersect.eps")
        fig.clear()

        #stat_frame = pd.DataFrame(
        #    {'avg_oall': avg_oall, 'num_cells': num_cells, 'var': var, 'varmean': varmean, 'rankmean': rankmean,
        #     'rankrat': rankrat, 'ranknum': ranknum})

        # plot mean vs variance to mean ratio

        fs = 14
        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_subplot(111)
        ax.set_title(label + clabel + glabel + "overlap with HVG2k set", fontsize=fs)
        ax.set_ylabel("variance to mean ratio")
        ax.set_xlabel("average expression")
        ax.loglog()
        sc1 = plt.scatter(stat_frame_all.avg_oall, stat_frame_all.varmean, c='black', s=80, edgecolors='None',
                          marker=".", label='all genes')
        sc1 = plt.scatter(stat_frame_gg.avg_oall, stat_frame_gg.varmean, c='deeppink', s=80, edgecolors='None',
                          marker=".", label=glabel)
        sc1 = plt.scatter(stat_frame_ddgs.avg_oall, stat_frame_ddgs.varmean, c='deepskyblue', s=80, edgecolors='None',
                          marker=".", label=title_list2[1])  # stat frame DDGs is ref group
        sc1 = plt.scatter(stat_frame_int.avg_oall, stat_frame_int.varmean, c='gold', s=80, edgecolors='None',
                          marker=".", label='both')
        ax.legend()
        plt.tight_layout()
        plt.savefig(large_root + "/" + label + glabel + clabel + "scatter_set_intersect_hvg2k_dispersion.png")
        plt.savefig(large_root + "/" + label + glabel + clabel + "scatter_set_intersect_hvg2k_dispersion.eps")
        fig.clear()
