import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.colors
import scipy.stats as ss
import seaborn as sns
from kneed import KneeLocator
from scipy.stats import variation


large_root = "/home/breanne/sourcedata/limb_bud/resample"
print("Reading data...")
dnum = 1 #number of density graphs to print


#e13 = pd.read_csv(large_root + r"/e13_limbbuddrop0.csv", index_col=0) #gene by cell
#e12 = pd.read_csv(large_root + r"/e12_limbbuddrop0.csv", index_col=0)
#e11 = pd.read_csv(large_root + r"/e11_limbbuddrop0.csv", index_col=0)
#e13_5 = pd.read_csv(large_root + r"/e13_5_limbbuddrop0.csv", index_col=0)
#e14 = pd.read_csv(large_root + r"/e14_limbbuddrop0.csv", index_col=0)
#e13b = pd.read_csv(large_root + r"/e13b_limbbuddrop0.csv", index_col=0)
#e15 = pd.read_csv(large_root + r"/e15_limbbuddrop0.csv", index_col=0)
#e10_5 = pd.read_csv(large_root + r"/e10_5_limbbuddrop0.csv", index_col=0)

label = 'e13_'
rep_list = [1,2,3,4,5,6,7,8,9,10]
fra_sig = []
num_genes = []
#num_counts = []
for rep in rep_list:
    rep = str(rep)
    pvals = pd.read_csv(large_root + "/pvals/rep" + rep + "_e13_5000_gene_params.tx.csv", index_col= 0, header = None, names = ["Pval"])

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

    #nonsig = bhv.index[bhv['sigp'] == 0].tolist()
    #DDG_list = data2.drop(nonsig, 'index') ############### needs to be data2
    #DDG_list.to_csv(large_root + "/" + label + "DDG_data.csv")

    num_sig_pvals = np.sum(bhv['sigp'])
    f_sig = num_sig_pvals/pvals.shape[0] #of genes that are expressed at all in data
    frac_sig = str(f_sig)
    print(frac_sig)
    num_g = str(pvals.shape[0])

    #num_counts.append(f_sig)
    fra_sig.append(f_sig)
    num_genes.append(num_g)
samp_data = pd.DataFrame({'fract_sig': fra_sig, 'num_genes': num_genes})
samp_data.to_csv(large_root + "/" + label + "resamp_data.csv")
