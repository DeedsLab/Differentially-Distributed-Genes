import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.colors
import scipy.stats as ss
import seaborn as sns
from kneed import KneeLocator
from scipy.stats import variation


large_root = "/home/breannesparta/ddgs/A20_3T3/cellgene/"
print("Reading data...")

label = 'A20_3T3_p90_rep'
#cell_list = [10000]
#rep_list = [1]
#total_cells = [80985]
rep_list = [1,2,3,4,5,6,7,8,9,10]
total_cells = [21123, 21123, 21123, 21123, 21123, 21123, 21123,21123, 21123, 21123, 21123, 21123]

#data2 = pd.read_csv("/media/Seagate_Repo/breanne/resample/zheng10k.csv", index_col=0) #must be transposed
#data2 = data2.T

#num_counts = []
for tcell in total_cells:
    fra_sig = []
    num_genes = []
    for rep in rep_list:
        rep = str(rep)
        ng = str(tcell)
        #data2 = pd.read_csv(sample_root + "/_rep_" + rep + "zheng9_resamp.csv", index_col=0)
        #data2 = data2.T
        pvals = pd.read_csv(large_root + "downsample/a203t3_p09_rep" + rep + "_pvals.txt", index_col= 0, header = None, names = ["Pval"])
        #pvals = pvals.loc[data2.index]
        print(pvals)
        print("files updated")

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
        sig_genes.to_csv(large_root + "/" + label + rep + "significant_genes.csv")

        nonsig = bhv.index[bhv['sigp'] == 0].tolist()
        #DDG_list = data2.drop(nonsig, 'index') ############### needs to be data2
        #DDG_list.to_csv("/media/timothyhamilton/Seagate Desktop Drive1/breanne/ddgs/zheng_resampled/ddgs/" + label + cell + "_rep_" + rep + "_DDGs.csv")

        num_sig_pvals = np.sum(bhv['sigp'])
        f_sig = num_sig_pvals/pvals.shape[0] #of genes that are expressed at all in data
        frac_sig = str(f_sig)
        print(frac_sig)
        num_g = str(pvals.shape[0])

        #num_counts.append(f_sig)
        fra_sig.append(f_sig)
        num_genes.append(num_g)
    samp_data = pd.DataFrame({'fract_sig': fra_sig, 'num_genes': num_genes})
    samp_data.to_csv(large_root + "/" + label + "_resamp_data.csv")
