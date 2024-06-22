import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.colors
import scipy.stats as ss
import seaborn as sns
from kneed import KneeLocator
from scipy.stats import variation
import scipy.stats as st


large_root = "/home/breannesparta/ddgs/zhengnine"
print("Reading data...")
dnum = 1 #number of density graphs to print

e10_5 = pd.read_csv(large_root + r"/zheng_reps_100_resamp_data.csv", index_col=0)
e11 = pd.read_csv(large_root + r"/zheng_reps_500_resamp_data.csv", index_col=0)
e12 = pd.read_csv(large_root + r"/zheng_reps_1000_resamp_data.csv", index_col=0)
e13 = pd.read_csv(large_root + r"/zheng_reps_2000_resamp_data.csv", index_col=0) #gene by cell
e13b = pd.read_csv(large_root + r"/zheng_reps_3000_resamp_data.csv", index_col=0)
e13_5 = pd.read_csv(large_root + r"/zheng_reps_4000_resamp_data.csv", index_col=0)
e14 = pd.read_csv(large_root + r"/zheng_reps_5000_resamp_data.csv", index_col=0)
e15 = pd.read_csv(large_root + r"/zheng_reps_6000_resamp_data.csv", index_col=0)
e16 = pd.read_csv(large_root + r"/zheng_reps_7000_resamp_data.csv", index_col=0) #gene by cell
e17 = pd.read_csv(large_root + r"/zheng_reps_8000_resamp_data.csv", index_col=0)
e18 = pd.read_csv(large_root + r"/zheng_reps_9000_resamp_data.csv", index_col=0)
e19 = pd.read_csv(large_root + r"/zheng_reps_10000_resamp_data.csv", index_col=0)

label = 'new_zheng10k'
data_list = [e10_5,e11,e12,e13,e13b,e13_5,e14,e15,e16,e17,e18,e19]
label_list = ['100','500','1000','2000','3000','4000','5000','6000', '7000', '8000', '9000', '10000']
fract_data = []
ngene_data = []
CI1 = []
CI2 = []
CI3 = []
CI4 = []
for data, cname in zip(data_list,label_list):
    fra_sig = data.fract_sig
    num_genes = data.num_genes
    fmean = np.mean(fra_sig)
    gmean = np.mean(num_genes)
    frac_interval = st.t.interval(alpha=0.95, df=len(fra_sig) - 1, loc=np.mean(fra_sig), scale=st.sem(fra_sig))
    gene_interval = st.t.interval(alpha=0.95, df=len(num_genes) - 1, loc=np.mean(num_genes), scale=st.sem(num_genes))
    CI1.append(fmean-frac_interval[0])
    CI2.append(frac_interval[1]-fmean)
    CI3.append(gmean-gene_interval[0])
    CI4.append(gene_interval[1]-gmean)
    fract_data.append(fmean)
    ngene_data.append(gmean)
samp_data = pd.DataFrame({'label_list': label_list, 'mean_sig': fract_data,'CIm1': CI1, 'CIm2': CI2, 'num_genes': ngene_data, 'CIg1': CI3, 'CIg2': CI4,})
samp_data.to_csv("/home/breannesparta/ddgs/zhengnine/" + label + "agg_resamp_data_10k.csv")

fs = 12
fig = plt.figure(figsize=(5, 8))
ax = fig.add_subplot(111)
ax = samp_data.plot.bar(x='label_list', y='mean_sig', yerr=[samp_data.CIm1, samp_data.CIm2], rot=0, width=0.75)
ax.set_title("Zheng: fraction significant DDGs per number of cells", fontsize=fs)
ax.set_ylabel("fraction ddgs")
ax.set_xlabel("n cells per type")
#ax.set_ylim(0, 1)
plt.grid()
plt.savefig("/home/breannesparta/ddgs/zhengnine/" + label + "frac_ddgs.png")
plt.savefig("/home/breannesparta/ddgs/zhengnine/" + label + "frac_ddgs.eps")

fs = 12
fig = plt.figure(figsize=(5, 8))
ax = fig.add_subplot(111)
ax = samp_data.plot.bar(x='label_list', y='num_genes', yerr=[samp_data.CIg1, samp_data.CIg2], rot=0, width=0.75)
ax.set_title("Zheng: gene count per sample", fontsize=fs)
ax.set_ylabel("num genes")
ax.set_xlabel("n cells per type")
#ax.set_ylim(0, 1)
plt.grid()
plt.savefig("/home/breannesparta/ddgs/zhengnine/" + label + "num_genes.png")
plt.savefig("/home/breannesparta/ddgs/zhengnine/" + label + "num_genes.eps")
