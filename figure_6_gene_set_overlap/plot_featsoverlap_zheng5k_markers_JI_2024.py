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

large_root = r"/home/breanne/ddgs"
label = 'zheng5k_plusmarkers'

total_data = pd.read_csv(r"/home/breanne/ddgs/zheng5k_dropzeros.csv", index_col=0, usecols =[0,1]) #gene x cells
markers = pd.read_csv(r"/home/breanne/ddgs/zhengnine/zheng_Panglao_markergenes.csv", index_col=1) #gene list
#cpm_data = pd.read_csv(r"/media/timothyhamilton/My Passport1/Tim_Hamilton/Zheng_FACS_pbmcs/zheng9_5k/zheng5k_cpmlog.csv", index_col=0) #cell x genes? use rows
#WC = pd.read_csv(r"/media/timothyhamilton/My Passport1/Tim_Hamilton/Zheng_FACS_pbmcs/zheng9_5k/zheng9_5k_WC.csv", index_col=0, usecols = [0,1]) #g x c
WC = pd.read_csv(r"/home/breanne/ddgs/zheng5k_oldWilcox_bf8k.csv", index_col=0, nrows=1) #c x g
HVG = pd.read_csv(r"/media/timothyhamilton/My Passport2/Tim_Hamilton/Zheng_FACS_pbmcs/zheng9_5k/zheng5k_HVGs.csv", index_col=0,skiprows=[1156],usecols = [0,1]) #remoe the extra indeces that come up from the duplicate gene
# 'Unnamed: 1156' - is index of 'NA-1' - gene, which is probably the duplicate total data gene that has same name when - are converted to .
DDG = pd.read_csv(r"/media/timothyhamilton/My Passport2/Tim_Hamilton/Zheng_FACS_pbmcs/zheng9_5k/zheng9_5k_DDGs.csv", index_col=0, usecols = [0,1])# g x c
townes = pd.read_csv(r"/media/timothyhamilton/My Passport2/Tim_Hamilton/Zheng_FACS_pbmcs/zheng9_5k/zheng5k_dropzeros_SCRY.csv",index_col = 0,nrows=1)
mdrop = pd.read_csv(r"/media/timothyhamilton/My Passport2/Tim_Hamilton/Zheng_FACS_pbmcs/zheng9_5k/zheng5k_dropzeros_M3Drop.csv",index_col = 0,nrows=1)
nbdrop = pd.read_csv(r"/media/timothyhamilton/My Passport2/Tim_Hamilton/Zheng_FACS_pbmcs/zheng9_5k/zheng5k_dropzeros_NBUMI.csv", index_col=0, nrows = 1)#c x g
nbdisp = pd.read_csv(r"/media/timothyhamilton/My Passport2/Tim_Hamilton/Zheng_FACS_pbmcs/zheng9_5k/zheng5k_dropzeros_NBUMI_High_Var.csv", index_col=0, nrows = 1) # c x g

marker_first = rename(markers.index.tolist())
WC_set= rename(WC.columns.tolist())
all_set= rename(total_data.index.tolist())
DDG_set= rename(DDG.index.tolist())
M3Drop_set= rename(mdrop.columns.tolist())
Townes_set= rename(townes.columns.tolist())
HVG_first = rename(HVG.index.tolist())
nbdrop_set= rename(nbdrop.columns.tolist())
nbdisp_set = rename(nbdisp.columns.tolist())

HVG_set = [x for x in HVG_first if x in all_set]
marker_set = [x for x in marker_first if x in all_set]
print(HVG_set)
print(marker_set)

title_list = ["all_genes","marker genes","WC","DDGs","M3Drop","Townes","NBdrop","NBdisp","HVGs"]
gene_list = [all_set, marker_set, WC_set, DDG_set, M3Drop_set, Townes_set,nbdrop_set, nbdisp_set,HVG_set]

#cell_headers = np.array(total_data.index)
#total_data = total_data.values
#total_data= pd.DataFrame(total_data, index = cell_headers, columns = gene_list[0])
title_list2 = title_list[1:]
gene_list2 = gene_list[1:]

for gene_group, glabel in zip(gene_list2, title_list2):
    print(glabel)
    print(len(gene_group))

nfeatures = [len(featureset) for featureset in gene_list2]
fig = plt.figure(figsize=(12, 8))
ax1 = fig.add_subplot(111)
nfeatures = pd.DataFrame({'gene_set': title_list2, 'AJI': nfeatures})
nfeatures.to_csv(large_root + "/" + label + "nfeatures.csv")
fs = 12
ax1 = nfeatures.plot.bar(x='gene_set', y='AJI', rot=0, width = 0.75)
ax1.set_ylabel("number of genes")
ax1.set_xlabel("feature set")
ax1.set_title("number of genes per set")
plt.grid()
plt.savefig(large_root + "/" + label + "n_featurse.eps")


#calc jaccard index with WC group
jindex = [WC_similarity(list1) for list1 in gene_list2] #marker set in this
fig = plt.figure(figsize=(12, 8))
ax1 = fig.add_subplot(111)
jindex_frame = pd.DataFrame({'gene_set': title_list2, 'AJI': jindex})
jindex_frame.to_csv(large_root + "/" + label + "jaccardindex_w_markergenes.csv")
fs = 12
ax1 = jindex_frame.plot.bar(x='gene_set', y='AJI', rot=0, width = 0.75)
ax1.set_ylabel("Jaccard Index")
ax1.set_xlabel("set compared to marker genes")
ax1.set_title("proportional overlap with marker genes")
ax1.set_ylim(0,1)
plt.grid()
plt.savefig(large_root + "/" + label + "jaccardI_w_markers.eps")
plt.savefig(large_root + "/" + label + "jaccardI_w_markers.png")

print("files updated")
