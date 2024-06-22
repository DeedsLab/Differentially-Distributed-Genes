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

large_root = r"/media/timothyhamilton/My Passport1/Tim_Hamilton/Zheng_FACS_pbmcs/zheng9_5k"
label = 'zheng5k'

total_data = pd.read_csv(r"/media/timothyhamilton/My Passport1/Tim_Hamilton/Zheng_FACS_pbmcs/zheng9_5k/zheng5k_dropzeros.csv", index_col=0) #gene x cells use cols
cpm_data = pd.read_csv(r"/media/timothyhamilton/My Passport1/Tim_Hamilton/Zheng_FACS_pbmcs/zheng9_5k/zheng5k_cpmlog.csv", index_col=0) #cell x genes? use rows
WC = pd.read_csv(r"/media/timothyhamilton/My Passport1/Tim_Hamilton/Zheng_FACS_pbmcs/zheng9_5k/zheng9_5k_WC.csv", index_col=0, usecols = [0,1]) #g x c
HVG = pd.read_csv(r"/media/timothyhamilton/My Passport1/Tim_Hamilton/Zheng_FACS_pbmcs/zheng9_5k/zheng5k_HVGs.csv", index_col=0,skiprows=[1156],usecols = [0,1]) #remoe the extra indeces that come up from the duplicate gene
# 'Unnamed: 1156' - is index of 'NA-1' - gene, which is probably the duplicate total data gene that has same name when - are converted to .
DDG = pd.read_csv(r"/media/timothyhamilton/My Passport1/Tim_Hamilton/Zheng_FACS_pbmcs/zheng9_5k/zheng9_5k_DDGs.csv", index_col=0, usecols = [0,1])# g x c
townes = pd.read_csv(r"/media/timothyhamilton/My Passport1/Tim_Hamilton/Zheng_FACS_pbmcs/zheng9_5k/zheng5k_dropzeros_SCRY.csv",index_col = 0,nrows=1)
mdrop = pd.read_csv(r"/media/timothyhamilton/My Passport1/Tim_Hamilton/Zheng_FACS_pbmcs/zheng9_5k/zheng5k_dropzeros_M3Drop.csv",index_col = 0,nrows=1)
nbdrop = pd.read_csv(r"/media/timothyhamilton/My Passport1/Tim_Hamilton/Zheng_FACS_pbmcs/zheng9_5k/zheng5k_dropzeros_NBUMI.csv", index_col=0, nrows = 1)#c x g
nbdisp = pd.read_csv(r"/media/timothyhamilton/My Passport1/Tim_Hamilton/Zheng_FACS_pbmcs/zheng9_5k/zheng5k_dropzeros_NBUMI_High_Var.csv", index_col=0, nrows = 1) # c x g
WC = WC.T
DDG = DDG.T
HVG = HVG.T
total_data = total_data.T

all_set= set(np.array(total_data.columns))
DDG_set= set(np.array(DDG.columns))
WC_set= set(np.array(WC.columns))
M3Drop_set= set(np.array(mdrop.columns))
Townes_set= set(np.array(townes.columns))
HVG_set = set(np.array(HVG.columns))
nbdrop_set= set(np.array(nbdrop.columns))
nbdisp_set = set(np.array(nbdisp.columns))

title_list = ["all_genes","WC","DDGs","HVGs","Townes","NBdrop"]
val_list = [all_set, WC_set, DDG_set, HVG_set, Townes_set,nbdrop_set]
gene_list = [rename(A) for A in val_list]

cell_headers = np.array(total_data.index)
total_data = total_data.values
total_data= pd.DataFrame(total_data, index = cell_headers, columns = gene_list[0])
title_list2 = title_list[1:]
gene_list2 = gene_list[1:]


#venn across all groups# max venn is 6
val_dict= dict(list(zip(title_list2,gene_list2)))
fig = plt.figure(figsize=(15, 8))
ax1 = fig.add_subplot(111)
ax1.set_title("Overlap of Gene_Sets")
venn(val_dict, cmap="viridis", fontsize=12, legend_loc="upper left", ax=ax1)
ax1.set_rasterized(True)
plt.savefig(large_root + "/" + label + "gene_overlap_NB2.eps")
plt.savefig(large_root + "/" + label + "gene_overlap_NB2.png")