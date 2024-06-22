import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.colors
import scipy.stats as ss
import seaborn as sns
from kneed import KneeLocator
from scipy.stats import variation
import scipy.stats as st


large_root = "/home/breannesparta/ddgs/zheng9"
print("Reading data...")
dnum = 1 #number of density graphs to print

p5data = pd.read_csv(large_root + r"/downsample/zheng9p05_newresamp_data_overlap_DDGp5.csv", index_col=0)
p9data = pd.read_csv(large_root + r"/downsample/zheng9p09_newresamp_data_overlap_DDGp5.csv", index_col=0)
allgenes = pd.read_csv(large_root + "/cellgene/zheng9pcap5DDG.csv", index_col=0, nrows=2)# c x g
print(allgenes.head())
allgenes = allgenes.T
#WC_names= WC.index.tolist()
#WC_set = rename(WC_names)
#print(len(WC_set))
#allgenes_set= rename(allgenes.index.tolist())
ngenes = len(allgenes.index.tolist())
print(ngenes)

label = 'zheng95k_downsample'
data_list = [p5data,p9data]
label_list = ['p50','p90']
fract_data = []
CI1 = []
CI2 = []
#CI3 = []
#CI4 = []
for data, cname in zip(data_list,label_list):
    fra_sig = (data.DDG_overlap)/ngenes
    print(fra_sig)
    #num_genes = data.WC_overlap
    fmean = np.mean(fra_sig)
    #gmean = np.mean(num_genes)
    frac_interval = st.t.interval(alpha=0.95, df=len(fra_sig) - 1, loc=np.mean(fra_sig), scale=st.sem(fra_sig))
    #gene_interval = st.t.interval(alpha=0.95, df=len(num_genes) - 1, loc=np.mean(num_genes), scale=st.sem(num_genes))
    CI1.append(fmean-frac_interval[0])
    CI2.append(frac_interval[1]-fmean)
    #CI3.append(gmean-gene_interval[0])
    #CI4.append(gene_interval[1]-gmean)
    fract_data.append(fmean)
    #ngene_data.append(gmean)
samp_data = pd.DataFrame({'label_list': label_list, 'frac_sig': fract_data,'CIm1': CI1, 'CIm2': CI2})
#samp_data = pd.DataFrame({'label_list': label_list, 'num_genes': ngene_data, 'CIg1': CI3, 'CIg2': CI4,})
samp_data.to_csv(large_root + "/downsample/" + label + "agg_resamp_data.csv")

fs = 12
fig = plt.figure(figsize=(5, 8))
ax = fig.add_subplot(111)
ax = samp_data.plot.bar(x='label_list', y='frac_sig', yerr=[samp_data.CIm1, samp_data.CIm2], rot=0, width=0.75)
ax.set_title("zheng9: DDG overlap in resampled data", fontsize=fs)
ax.set_ylabel("fraction ddg recovered in downsampled data")
ax.set_xlabel("capture probability")
ax.set_ylim(0, 1)
plt.grid()
#plt.savefig(large_root + "/" + label + "ddgs_over.png")
plt.savefig(large_root + "/downsample/" + label + "ddgs_over.eps")


"""fs = 12
fig = plt.figure(figsize=(5, 8))
ax = fig.add_subplot(111)
ax = samp_data.plot.bar(x='label_list', y='num_genes', yerr=[samp_data.CIg1, samp_data.CIg2], rot=0, width=0.75)
ax.set_title("Zheng: DDG overlap with WC in resampled data", fontsize=fs)
ax.set_ylabel("fraction ddg-WC overlap")
ax.set_xlabel("n cells per type")
#ax.set_ylim(0, 1)
plt.grid()
plt.savefig(large_root + "/" + label + "WC_over.png")
plt.savefig(large_root + "/" + label + "WC_over_.eps")"""
