import pandas as pd
import numpy as np
import scanpy as sc
from sklearn.metrics.cluster import adjusted_rand_score
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.neighbors import NearestNeighbors
from sklearn.neighbors import DistanceMetric
from sklearn.decomposition import PCA
from kneed import KneeLocator as kl


def rename(A):
    # for gene list A return gene names with spaces instead of - or . characters
    newlist = []
    for gene in A:
        newgene = str(gene)
        #a = newgene.replace('.', '').replace('-', '') #zheng has a gene where if you remove both . and -, two genes have the same name
        # probably where nan comes from in hvg data
        a = newgene.replace('.', '-') #R changes - to .
        newlist.append(a)
    all_headers = np.array(newlist)
    return (all_headers)

large_root = "/home/breanne/sourcedata/Zheng9"
label = 'zheng5k_mergTls'

#usec = list(range(0,101))
total_data = pd.read_csv(r"/home/breanne/sourcedata/Zheng9/zheng5k_dropzeros.csv", index_col=0) #gene x cells use cols
cpm_data = pd.read_csv(r"/home/breanne/sourcedata/Zheng9/zheng5k_cpmlog.csv", index_col=0) #cell x genes? use rows
WC = pd.read_csv(r"/home/breanne/sourcedata/Zheng9/zheng9_5k_WC.csv", index_col=0, usecols = [0,1]) #g x c
#HVG = pd.read_csv(r"/home/breanne/sourcedata/Zheng9/zheng5k_HVGs.csv", index_col=0,skiprows=[1156],usecols = [0,1]) #remoe the extra indeces that come up from the duplicate gene
## 'Unnamed: 1156' - is index of 'NA-1' - gene, which is probably the duplicate total data gene that has same name when - are converted to .
DDG = pd.read_csv(r"/home/breanne/sourcedata/Zheng9/zheng9_5k_DDGs.csv", index_col=0, usecols = [0,1])# g x c
#townes = pd.read_csv(r"/home/breanne/sourcedata/Zheng9/zheng5k_dropzeros_SCRY.csv",index_col = 0,nrows=1)
#mdrop = pd.read_csv("/media/timothyhamilton/data1/Tim_Hamilton/Zheng_FACS_pbmcs/zheng9_5k/zheng5k_dropzeros_M3Drop.csv",index_col = 0,nrows=1)
#nbdrop = pd.read_csv("/media/timothyhamilton/data1/Tim_Hamilton/Zheng_FACS_pbmcs/zheng9_5k/zheng5k_dropzeros_NBUMI.csv", index_col=0, nrows = 1)#c x g
#nbdisp = pd.read_csv("/media/timothyhamilton/data1/Tim_Hamilton/Zheng_FACS_pbmcs/zheng9_5k/zheng5k_dropzeros_NBUMI_High_Var.csv", index_col=0, nrows = 1) # c x g
WC = WC.T
DDG = DDG.T
total_data = total_data.T
#HVG = HVG.T

#clean gene names so they can be used to subset on all data
WC_set= rename(WC.columns.tolist())
all_set= rename(total_data.columns.tolist())
DDG_set= rename(DDG.columns.tolist())
#M3Drop_set= rename(mdrop.columns.tolist())
#Townes_set= rename(townes.columns.tolist())
#HVG_first = rename(HVG.columns.tolist())
#nbdrop_set= rename(nbdrop.columns.tolist())
#nbdisp_set = rename(nbdisp.columns.tolist())
#HVG_set = [x for x in HVG_first if x in all_set] #strange NA-1 space maybe, that doesnt exist in all data

#title_list = ["WC","DDGs","M3Drop","Townes","NBdrop","NBdisp","HVGs","all_genes"]
#gene_list = [WC_set, DDG_set, M3Drop_set, Townes_set, nbdrop_set, nbdisp_set, HVG_set, all_set]

#title_list = ["WC","DDGs","HVGs","all_genes"]
#gene_list = [WC_set, DDG_set, HVG_set, all_set]

#HVGs are giving value error WHY idk
title_list = ["WC","DDGs","all_genes"]
gene_list = [WC_set, DDG_set,all_set]

#group confusion groups into the same index
facs_headers = np.array(total_data.index)  # headers h is index, so data is cell x gene
top_clusterin = pd.DataFrame(index=facs_headers)
top_clusterin['FACS'] = facs_headers
top_clusterin = top_clusterin.replace('^.*b_cells.*$', '1', regex=True)
top_clusterin = top_clusterin.replace('^.*cd4_t_helper.*$', '2', regex=True) #
top_clusterin = top_clusterin.replace('^.*cd56_nk.*$', '3', regex=True)
top_clusterin = top_clusterin.replace('^.*cytotoxic_t.*$', '4', regex=True) ##
top_clusterin = top_clusterin.replace('^.*memory_t.*$', '4', regex=True) ##
top_clusterin = top_clusterin.replace('^.*monocytes.*$', '6', regex=True)
top_clusterin = top_clusterin.replace('^.*naive_cytotoxic.*$', '4', regex=True) ##
top_clusterin = top_clusterin.replace('^.*naive_t.*$', '2', regex=True) #
top_clusterin = top_clusterin.replace('^.*regulatory_t.*$', '2', regex=True) #
PC_clusterin = pd.DataFrame(index=facs_headers)
PC_clusterin['FACS'] = top_clusterin.FACS

total_data = total_data.values
total_data= pd.DataFrame(total_data, index = facs_headers, columns = all_set)
cpm_data = cpm_data.values
cpm_data= pd.DataFrame(cpm_data, index = facs_headers, columns = all_set)
#eq_data = total_data.values
#eq_data = equalizer(eq_data)
#eq_data = pd.DataFrame(eq_data, index = facs_headers, columns = gene_list[0])

data_list = [total_data, cpm_data]
countname = ['_raw_counts_', '_cpm_log_counts_']

#cycle through different normalizations
for data_group, cname in zip(data_list,countname):

    # cycle through different feature sets
    for glabel, gname in zip(gene_list,title_list):

        #this is identitical to plot feats overlap that works with HVGs
        test_group = pd.Series(list(glabel))
        new_data = data_group.loc[:, test_group]
        new_data = new_data.fillna(0) #there shouldnt be any though

        data = sc.AnnData(new_data)
        data.obs['facs'] = top_clusterin.FACS
        sc.tl.pca(data, n_comps=50)
        sc.tl.tsne(data, n_pcs=50)

        # make eps files I want with colors by avoiding scanpy
        tsne_cords = pd.DataFrame(data.obsm["X_tsne"], index=facs_headers)
        tsne_cords.rename(columns={tsne_cords.columns[0]: "tsne1"}, inplace=True)
        tsne_cords.rename(columns={tsne_cords.columns[1]: "tsne2"}, inplace=True)
        tsne_cords['FACS'] = top_clusterin['FACS']

        # there are 9 LV, 5 FACS ############ will have to edit this
        LV_colors = ['dodgerblue', 'darkgreen', 'firebrick', 'rebeccapurple', \
                        'cyan', 'grey', 'darkorange', 'magenta']
        FACS_colors = ['dodgerblue', 'darkgreen', 'firebrick', 'rebeccapurple', \
                        'cyan', 'grey', 'darkorange', 'magenta']
        fs = 10
        fig = plt.figure(figsize=(5, 8))
        ax = fig.add_subplot(111)
        sns.set_palette(sns.color_palette(FACS_colors))
        sns.lmplot(x="tsne1", y="tsne2", data=tsne_cords, fit_reg=False, hue='FACS', legend=None, scatter_kws={"s": 2})
        ax = plt.gca()
        ax.set_title(label + gname + cname + " colored byFACS ", fontsize=fs)
        ax.set_ylabel("tSNE_2")
        ax.set_xlabel("tSNE_1")
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.tight_layout()
        plt.savefig(large_root + "/" + label + gname + cname + "FACS_TSNE2.png")
        plt.savefig(large_root + "/" + label + gname + cname + "FACS_TSNE2.eps")

