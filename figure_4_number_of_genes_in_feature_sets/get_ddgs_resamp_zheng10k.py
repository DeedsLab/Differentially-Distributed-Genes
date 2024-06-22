import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.colors
import scipy.stats as ss
import seaborn as sns
from kneed import KneeLocator
from scipy.stats import variation


large_root = "/home/breanne/sourcedata/Zheng9/resample"
print("Reading data...")

label = 'zheng_'
cell_list = [10000]
rep_list = [1]
total_cells = [80985]
#cell_list = [100,500,1000,2000,3000,4000,5000,6000, 7000, 8000, 9000, 10000]
#rep_list = [1,2,3,4,5,6,7,8,9,10]
#total_cells = [809, 4049, 8098, 16197, 24295, 32394, 40492,48591, 56689, 64788, 72886, 80985]
label = 'zheng10k'
data2 = pd.read_csv("/media/Seagate_Repo/breanne/resample/zheng10k.csv", index_col=0) #must be transposed
#data2 = data2.T

fra_sig = []
num_genes = []
#num_counts = []
for cell, tcell in zip(cell_list, total_cells):
    for rep in rep_list:
        rep = str(rep)
        cell = str(cell)
        ng = str(tcell)
        pvals = pd.read_csv(large_root + "/new_resamp10k_pvals/ep_" + rep + "zheng_resamp" + ng + "gene_params.tx.csv", header = None, names = ["gene", "Pval"])
        #pvals = pd.read_csv(large_root + "/resample_pvals/" + cell + "_rep_" + rep + "zheng_resamp" + ng + "gene_params.tx.csv", index_col= 0, header = None, names = ["Pval"])

        dropzero = data2[np.all(data2 == 0, axis=1)].index
        data2 = data2.drop(dropzero, 'index')  # drop genes that are not detected in any cells
        pvals = pvals.drop(dropzero, 'index')  # drop genes that are not detected in any cells

        gene_headers = np.array(data2.index)
        pvals = pvals.drop(columns=['gene'])
        pvals.index = gene_headers

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
        sig_genes = pd.DataFrame(sgnames)
        sig_genes.to_csv(large_root + "/" + label + cell + rep + "significant_genes.csv")

        nonsig = bhv.index[bhv['sigp'] == 0].tolist()
        DDG_list = data2.drop(nonsig, 'index') ############### needs to be data2
        DDG_list.to_csv("/media/Seagate_Repo/breanne/resample/" + label + "DDG_data.csv")

        num_sig_pvals = np.sum(bhv['sigp'])
        f_sig = num_sig_pvals/pvals.shape[0] #of genes that are expressed at all in data
        frac_sig = str(f_sig)
        print(frac_sig)
        num_g = str(pvals.shape[0])

        #num_counts.append(f_sig)
        fra_sig.append(f_sig)
        num_genes.append(num_g)
    #samp_data = pd.DataFrame({'fract_sig': fra_sig, 'num_genes': num_genes})
    #samp_data.to_csv(large_root + "/" + label + cell + "_resamp_data.csv")
