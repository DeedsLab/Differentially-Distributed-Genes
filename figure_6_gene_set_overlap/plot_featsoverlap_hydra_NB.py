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
total_data = pd.read_csv("/media/timothyhamilton/data1/Tim_Hamilton/Hydra/Hydra_full_Raw.csv", index_col=0, nrows =1) #cell x genes
#cpm_data = pd.read_csv("/media/timothyhamilton/data1/Tim_Hamilton/Hydra/Hydra_cpmlog.csv", index_col=0) #cell x genes
HVG = pd.read_csv("/media/timothyhamilton/data1/Tim_Hamilton/Hydra/Hydra_HVGs.csv", index_col=0, nrows=1) #c x g
DDG = pd.read_csv("/media/timothyhamilton/data1/Tim_Hamilton/Hydra/Hydra_DDG_data.csv", index_col=0, usecols = [0,1])# g x c
townes = pd.read_csv("/media/timothyhamilton/data1/Tim_Hamilton/Hydra/Hydra_Raw_Cell_Types_Townes.csv", index_col=0, nrows=1)#c x g
mdrop = pd.read_csv("/media/timothyhamilton/data1/Tim_Hamilton/Hydra/Hydra_M3Drop_Genes.csv", index_col=0, nrows=1)
nbdrop = pd.read_csv("/media/timothyhamilton/data1/Tim_Hamilton/R01_Ref_Datasets/Hydra_Results/Hydra_NBUMIDrop_Genes.csv", index_col=0, nrows=1)
nbdisp = pd.read_csv("/media/timothyhamilton/data1/Tim_Hamilton/R01_Ref_Datasets/Hydra_Results/Hydra_NBUMIDrop_Genes_High_Var.csv", index_col=0, nrows=1)

dropzero = total_data[np.all(total_data == 0, axis=1)].index
total_data = total_data.drop(dropzero, 'index') #drop genes that are not detected in any cells
#cpm_data = cpm_data.drop(dropzero, 'index') #drop genes that are not detected in any cells

nbdrop = nbdrop.drop(['Cell_Type'], axis = 1)
nbdisp = nbdisp.drop(['Cell_Type'], axis = 1)
townes = townes.drop(['Cell_Type'], axis = 1)
mdrop = mdrop.drop(['Cell_Type'], axis = 1)
DDG = DDG.T

all_set= set(np.array(total_data.columns))
DDG_set= set(np.array(DDG.columns))
M3Drop_set= set(np.array(mdrop.columns))
Townes_set= set(np.array(townes.columns))
HVG_set = set(np.array(HVG.columns))
nbdrop_set= set(np.array(nbdrop.columns))
nbdisp_set = set(np.array(nbdisp.columns))

title_list = ["all_genes","DDGs","HVGs","Townes","NBdrop"]
val_list = [all_set, DDG_set, HVG_set, Townes_set,nbdrop_set]
gene_list = [rename(S) for S in val_list]

cell_headers = np.array(total_data.index)
total_data = total_data.values
total_data= pd.DataFrame(total_data, index = cell_headers, columns = gene_list[0])
#cpm_data = cpm_data.values
#cpm_data= pd.DataFrame(cpm_data, index = cell_headers, columns = gene_list[0])
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
#title_list3 = title_list2[1:]
#gene_list3 = gene_list2[1:]
val_dict= dict(list(zip(title_list2,gene_list2)))
fig = plt.figure(figsize=(15, 8))
ax1 = fig.add_subplot(111)
ax1.set_title("Overlap of Gene_Sets")
venn(val_dict, cmap="viridis", fontsize=12, legend_loc="upper left", ax=ax1)
ax1.set_rasterized(True)
plt.savefig(large_root + "/" + label + "gene_overlap_NB.eps")
plt.savefig(large_root + "/" + label + "gene_overlap_NB.png")
