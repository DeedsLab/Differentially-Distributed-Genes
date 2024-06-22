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
    headers = np.array(newlist)
    return (headers)

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

large_root = "/media/timothyhamilton/data1/Tim_Hamilton/Hydra"
label = 'Hydra'
#AJD calc on cell by gene data
total_data = pd.read_csv("/media/timothyhamilton/data1/Tim_Hamilton/Hydra/Hydra_full_Raw.csv", index_col=0) #cell x genes
cpm_data = pd.read_csv("/media/timothyhamilton/data1/Tim_Hamilton/Hydra/Hydra_cpmlog.csv", index_col=0) #cell x genes
HVG = pd.read_csv("/media/timothyhamilton/data1/Tim_Hamilton/Hydra/Hydra_HVGs.csv", index_col=0, nrows=1) #c x g
DDG = pd.read_csv("/media/timothyhamilton/data1/Tim_Hamilton/Hydra/Hydra_DDG_data.csv", index_col=0, usecols = [0,1])# g x c
townes = pd.read_csv("/media/timothyhamilton/data1/Tim_Hamilton/Hydra/Hydra_Raw_Cell_Types_Townes.csv", index_col=0, nrows=1)#c x g
mdrop = pd.read_csv("/media/timothyhamilton/data1/Tim_Hamilton/Hydra/Hydra_M3Drop_Genes.csv", index_col=0, nrows=1)
nbdrop = pd.read_csv("/media/timothyhamilton/data1/Tim_Hamilton/R01_Ref_Datasets/Hydra_Results/Hydra_NBUMIDrop_Genes.csv", index_col=0, nrows=1)
nbdisp = pd.read_csv("/media/timothyhamilton/data1/Tim_Hamilton/R01_Ref_Datasets/Hydra_Results/Hydra_NBUMIDrop_Genes_High_Var.csv", index_col=0, nrows=1)

dropzero = total_data[np.all(total_data == 0, axis=1)].index
total_data = total_data.drop(dropzero, 'index') #drop genes that are not detected in any cells
cpm_data = cpm_data.drop(dropzero, 'index') #drop genes that are not detected in any cells

nbdrop = nbdrop.drop(['Cell_Type'], axis = 1)
nbdisp = nbdisp.drop(['Cell_Type'], axis = 1)
townes = townes.drop(['Cell_Type'], axis = 1)
mdrop = mdrop.drop(['Cell_Type'], axis = 1)
DDG = DDG.T

all_set= rename(total_data.columns.tolist())
DDG_set= rename(DDG.columns.tolist())
M3Drop_set= rename(mdrop.columns.tolist())
Townes_set= rename(townes.columns.tolist())
HVG_set = rename(HVG.columns.tolist())
nbdrop_set= rename(nbdrop.columns.tolist())
nbdisp_set = rename(nbdisp.columns.tolist())

title_list = ["all_genes","DDGs","HVGs","M3Drop","Townes","NBdrop","NBdisp"]
gene_list = [all_set, DDG_set, HVG_set, M3Drop_set, Townes_set,nbdrop_set,nbdisp_set]

cell_headers = np.array(total_data.index)
total_data = total_data.values
total_data= pd.DataFrame(total_data, index = cell_headers, columns = all_set)
cpm_data = cpm_data.values
cpm_data= pd.DataFrame(cpm_data, index = cell_headers, columns = all_set)
title_list2 = title_list[1:]
gene_list2 = gene_list[1:]

#calc jaccard index with WC group
jindex = [WC_similarity(list1) for list1 in gene_list2]
fig = plt.figure(figsize=(12, 8))
ax1 = fig.add_subplot(111)
jindex_frame = pd.DataFrame({'gene_set': title_list2, 'AJI': jindex})
jindex_frame.to_csv(large_root + "/" + label + "jaccardindex_w_DDGs.csv")
fs = 12
ax1 = jindex_frame.plot.bar(x='gene_set', y='AJI', rot=0, width = 0.75)
ax1.set_ylabel("Jaccard Index")
ax1.set_xlabel("set compared to DDGs")
ax1.set_title("proportional overlap with DDGs")
ax1.set_ylim(0,1)
plt.grid()
plt.savefig(large_root + "/" + label + "jaccardI_w_DDGs.eps")
plt.savefig(large_root + "/" + label + "jaccardI_w_DDGs.png")

#venn across all groups# max venn is 6
#needs to be a set, we already made this
"""val_dict= dict(list(zip(title_list2,gene_list2)))
fig = plt.figure(figsize=(15, 8))
ax1 = fig.add_subplot(111)
ax1.set_title("Overlap of Gene_Sets")
venn(val_dict, cmap="viridis", fontsize=12, legend_loc="upper left", ax=ax1)
ax1.set_rasterized(True)
plt.savefig(large_root + "/" + label + "gene_overlap.eps")
plt.savefig(large_root + "/" + label + "gene_overlap.png")"""

print("files updated")

data_list = [total_data, cpm_data]
cellname = ['_raw_counts_','_cpm_log_counts_']

for cell_group, clabel in zip(data_list,cellname):
    data2 = cell_group.T #gene x cell
    test_data = data2
    stat_frame_all = get_gene_stats(test_data)

    # get expectation for bionomial model
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
        # DDG set overlaps
        fig = plt.figure(figsize=(15, 8))
        ax = fig.add_subplot(111)
        venn2([set(gene_list[1]), set(gene_group)], set_labels=('DDGs', glabel))
        ax.set_rasterized(True)
        plt.savefig(large_root + "/" + label + clabel + "gene_overlap_DDG" + glabel + ".png")
        plt.savefig(large_root + "/" + label + clabel + "gene_overlap_DDG" + glabel + ".eps")

        #scatter gene sets with intersections
        test_group = pd.Series(list(gene_group))
        stat_frame_gg = stat_frame_all.loc[test_group, :]

        #intersect w DDGs
        ddg_group = pd.Series(list(gene_list[1]))
        stat_frame_ddgs = stat_frame_all.loc[ddg_group,:]

        set_2 = frozenset(gene_list[1])
        intersectWC = [x for x in gene_group if x in set_2]
        #intersectWC = gene_list[1].intersection(gene_group)
        #sets dont keep the order and learn about dif b/w view and copy
        i = len(intersectWC)
        print("intersect"+str(i))
        intersect_group = pd.Series(list(intersectWC))
        #stat_frame_int = stat_frame_all.filter(intersectWC, axis=0) #does filter not work?
        stat_frame_int = stat_frame_all.loc[intersect_group,:] #does filter not work?

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
        ax.set_title(label + clabel + glabel + "overlap with DDG set", fontsize=fs)
        ax.set_ylabel("number of cells")
        ax.set_xlabel("average mRNA count over all cells")
        ax.loglog()
        sc1 = plt.scatter(stat_frame_all.avg_oall, stat_frame_all.num_cells, c='black', s=80, edgecolors='None',
                          marker=".", label = 'all genes')
        sc1 = plt.scatter(stat_frame_gg.avg_oall, stat_frame_gg.num_cells, c='deeppink', s=80, edgecolors='None',
                          marker=".", label = glabel)
        sc1 = plt.scatter(stat_frame_ddgs.avg_oall, stat_frame_ddgs.num_cells, c='deepskyblue', s=80, edgecolors='None',
                          marker=".", label= 'DDGs')
        sc1 = plt.scatter(stat_frame_int.avg_oall, stat_frame_int.num_cells, c='gold', s=80, edgecolors='None',
                          marker=".", label='both')
        ax.legend()
        sc1 = plt.plot(ma, Yexp, linewidth=3)
        plt.tight_layout()
        plt.savefig(large_root + "/" + label + glabel + clabel+ "scatter_set_intersect_DDG.png")
        plt.savefig(large_root + "/" + label + glabel + clabel+ "scatter_set_intersect_DDG.eps")
        fig.clear()

        fs = 14
        fig = plt.figure(figsize=(18, 6))
        ax = fig.add_subplot(131)
        ax.set_title(label + clabel + glabel + "overlap with DDG set", fontsize=fs)
        ax.set_ylabel("average mRNA count over all cells")
        ax.set_xlabel("rank average count")
        sc1 = plt.scatter(stat_frame_all.rankmean, stat_frame_all.avg_oall, c='black', s=80, edgecolors='None',
                          marker=".", label='all genes')
        sc1 = plt.scatter(stat_frame_gg.rankmean, stat_frame_gg.avg_oall, c='deeppink', s=80, edgecolors='None',
                          marker=".", label=glabel)
        sc1 = plt.scatter(stat_frame_ddgs.rankmean, stat_frame_ddgs.avg_oall, c='deepskyblue', s=80, edgecolors='None',
                          marker=".", label='DDGs')
        sc1 = plt.scatter(stat_frame_int.rankmean, stat_frame_int.avg_oall, c='gold', s=80, edgecolors='None',
                          marker=".", label='both')
        ax.legend()
        plt.tight_layout()

        ax = fig.add_subplot(132)
        ax.set_title(label + clabel + glabel + "overlap with DDG set", fontsize=fs)
        ax.set_ylabel("number of cells mRNA is detected in")
        ax.set_xlabel("rank num cells")
        sc1 = plt.scatter(stat_frame_all.ranknum, stat_frame_all.num_cells, c='black', s=80, edgecolors='None',
                          marker=".", label='all genes')
        sc1 = plt.scatter(stat_frame_gg.ranknum, stat_frame_gg.num_cells, c='deeppink', s=80, edgecolors='None',
                          marker=".", label=glabel)
        sc1 = plt.scatter(stat_frame_ddgs.ranknum, stat_frame_ddgs.num_cells, c='deepskyblue', s=80, edgecolors='None',
                          marker=".", label='DDGs')
        sc1 = plt.scatter(stat_frame_int.ranknum, stat_frame_int.num_cells, c='gold', s=80, edgecolors='None',
                          marker=".", label='both')
        ax.legend()
        plt.tight_layout()

        ax = fig.add_subplot(133)
        ax.set_title(label + clabel + glabel + "overlap with DDG set", fontsize=fs)
        ax.set_ylabel("rank number of cells")
        ax.set_xlabel("rank average mRNA count over all cells")
        #ax.loglog()
        sc1 = plt.scatter(stat_frame_all.rankmean, stat_frame_all.ranknum, c='black', s=80, edgecolors='None',
                          marker=".", label='all genes')
        sc1 = plt.scatter(stat_frame_gg.rankmean, stat_frame_gg.ranknum, c='deeppink', s=80, edgecolors='None',
                          marker=".", label=glabel)
        sc1 = plt.scatter(stat_frame_ddgs.rankmean, stat_frame_ddgs.ranknum, c='deepskyblue', s=80, edgecolors='None',
                          marker=".", label='DDGs')
        sc1 = plt.scatter(stat_frame_int.rankmean, stat_frame_int.ranknum, c='gold', s=80, edgecolors='None',
                          marker=".", label='both')
        ax.legend()
        plt.tight_layout()
        plt.savefig(large_root + "/" + label + glabel + clabel + "scatter_set_intersect_DDG_rankmean_rcount.png")
        plt.savefig(large_root + "/" + label + glabel + clabel + "scatter_set_intersect_DDG_rankmean_rcount.eps")
        fig.clear()

        fs = 14
        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_subplot(111)
        ax.set_title(label + clabel + glabel + "overlap with DDG set", fontsize=fs)
        ax.set_ylabel("variance to mean ratio")
        ax.set_xlabel("average expression")
        ax.loglog()
        sc1 = plt.scatter(stat_frame_all.avg_oall, stat_frame_all.varmean, c='black', s=80, edgecolors='None',
        marker = ".", label = 'all genes')
        sc1 = plt.scatter(stat_frame_gg.avg_oall, stat_frame_gg.varmean, c='deeppink', s=80, edgecolors='None',
        marker = ".", label = glabel)
        sc1 = plt.scatter(stat_frame_ddgs.avg_oall, stat_frame_ddgs.varmean, c='deepskyblue', s=80, edgecolors='None',
        marker = ".", label = 'DDGs')
        sc1 = plt.scatter(stat_frame_int.avg_oall, stat_frame_int.varmean, c='gold', s=80, edgecolors='None',
        marker = ".", label = 'both')
        ax.legend()
        plt.tight_layout()
        plt.savefig(large_root + "/" + label + glabel + clabel + "scatter_set_intersect_DDG_dispersion.png")
        plt.savefig(large_root + "/" + label + glabel + clabel + "scatter_set_intersect_DDG_dispersion.eps")
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
    avg_oall = np.nanmean(test_data, 1)  # avg mRNA exp across all cells
    test_data[test_data == 0] = np.nan
    num_cells = test_data.count(axis=1)  # number of cells each mRNA is detected in
    perc_cells = [(x / NT) * 100 for x in num_cells]

    count_frame = pd.DataFrame({'avg_count': avg_oall, 'perc_cells': perc_cells})
    count_data = count_frame.melt().assign(x=glabel)

    ax = fig.add_subplot(1,len(title_list),c)
    ax.set_title(label + glabel , fontsize=fs)
    ax.set_ylabel("geneset")
    ax.set_xlabel("value")
    plt.ylim((-20, 100))
    ax.set(yscale="linear")
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
ax.set(yscale="log")
ax = sns.violinplot(data=avg_list, scale = "count")
ax.set(xticklabels=title_list)

ax = fig.add_subplot(122)
ax.set_title("percent cells mRNA is detected in for:" + label + "set", fontsize=fs)
ax.set_ylabel("percent cell count per mRNA")
ax.set_xlabel("geneset")
ax.set(yscale="log")
ax = sns.violinplot(data=perc_list, scale = "count")
ax.set(xticklabels=title_list)
plt.savefig(large_root + "/" + label + "compare_count_stats_nl.png")
plt.savefig(large_root + "/" + label + "compare_count_stats_nl.eps")