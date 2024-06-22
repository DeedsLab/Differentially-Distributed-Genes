import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.colors
import scipy.stats as ss
import seaborn as sns
from kneed import KneeLocator
from scipy.stats import variation
import scipy.stats as st

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

large_root = "/home/breannesparta/ddgs/zhengnine"
print("Reading data...")
WC = pd.read_csv("/home/breannesparta/ddgs/zheng5k_oldWilcox_bf8k.csv", index_col=0, nrows=1)
WC = WC.T
print(WC.head()) #genex x cells
#why was index col 1?
DDG = pd.read_csv(r"/media/timothyhamilton/Seagate Desktop Drive1/breanne/ddgs/zheng_resampled/ddgs/zheng_reps_10000_rep_1_DDGs.csv", index_col=0, usecols=[0])# g x c
WC_names= WC.index.tolist()
WC_set = rename(WC_names)
print(len(WC_set))
DDG_set= rename(DDG.index.tolist())
print(len(DDG_set))

label = 'zheng_5k'
cell_list = [100,500,1000,2000,3000,4000,5000,6000, 7000, 8000, 9000, 10000]
rep_list = [1,2,3,4,5,6,7,8,9,10]
total_cells = [809, 4049, 8098, 16197, 24295, 32394, 40492,48591, 56689, 64788, 72886, 80985]

#num_counts = []
for cell, tcell in zip(cell_list, total_cells):
    DDG_overlap = []
    WC_overlap = []
    for rep in rep_list:
        rep = str(rep)
        cell = str(cell)
        ng = str(tcell)
        pvals = pd.read_csv(large_root + "/pvals/" + cell + "_rep_" + rep + "zheng_resam_pvals.csv",
                            index_col=0, header=None, names=["Pval"])
        #10000_rep_10zheng_resam_pvals.csv
        #BH correction
        pcut = 0.01
        rankps = ss.rankdata(pvals, method='min')
        bhv = pd.DataFrame((rankps/pvals.shape[0])*pcut) ### is this number of genes
        bhv.columns = ['iQ/m']
        bhv.index =pvals.index
        bhv['rank'] = rankps
        bhv['pval'] = pvals
        #find largest pval that is smaller than iQm
        bhv['issig'] = np.where(bhv['pval']<bhv['iQ/m'],1,0)
        bhv['lowv'] = np.where(bhv['issig'] == 1,bhv['pval'],0)
        maxsig = np.max(bhv['lowv'])
        bhv['sigp'] = np.where(bhv['pval'] < maxsig,1,0)
        sgnames = bhv.index[bhv['sigp'] == 1].tolist()
        resamp_set = rename(sgnames)

        # scatter gene sets with intersections
        test_group = pd.Series(list(resamp_set))
        set_1 = frozenset(DDG_set.tolist())
        set_2 = frozenset(WC_set.tolist())
        print(set_2)

        h = len(set_2)
        print("WC genes" + str(h))

        intersectDDG = [x for x in test_group if x in set_1]
        intersectWC = [x for x in test_group if x in set_2]
        i = len(intersectDDG)
        #i = len(test_group)
        j = len(intersectWC)
        print("intersectDDGs" + str(i))
        print("intersectWC" + str(j))

        #num_counts.append(f_sig)
        DDG_overlap.append(i)
        WC_overlap.append(j)
    samp_data = pd.DataFrame({'DDG_overlap': DDG_overlap, 'WC_overlap': WC_overlap})
    #samp_data = pd.DataFrame({'WC_overlap': WC_overlap})
    samp_data.to_csv(large_root + "/" + label + cell + "_newresamp_data_overlap_WC_5kbf8kg.csv")
