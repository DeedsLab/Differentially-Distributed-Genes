import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.colors
import scipy.stats as ss
import seaborn as sns
from kneed import KneeLocator
from scipy.stats import variation


large_root = "/home/breanne/sourcedata/Zheng9"
print("Reading data...")
dnum = 1 #number of density graphs to print

label = 'zheng9_5k'
data = pd.read_csv("/home/breanne/sourcedata/Zheng9/zheng9_5k.csv") #must be transposed???
pvals = pd.read_csv("/home/breanne/sourcedata/Zheng9/zheng9_5k_Pvals.csv", header = None, names = ["gene", "Pval"])
data = data.set_index('Unnamed: 0')
gene_headers = np.array(data.index)
#data = data.T #?? confusion
pvals = pvals.drop(columns=['gene'])
pvals.index = gene_headers

cellgroup1 = data.filter(regex='^b',axis=1) #b_cells
cellgroup2 = data.filter(regex='^naive_cy',axis=1) #cytotoxic
cellgroup3 = data.filter(regex='^mono',axis=1) #monocytes
cellgroup4 = data.filter(regex='^regulatory',axis=1) #t
cellgroup5 = data.filter(regex='^cd4',axis=1) #thelper
cellgroup6 = data.filter(regex='^cd56',axis=1) #nk
cellgroup7 = data.filter(regex='^mem',axis=1) #memory t
cellgroup8 = data.filter(regex='^naive_t',axis=1) #naive t
cellgroup9 = data.filter(regex='^cytotoxic',axis=1) #cytotoxic_t
print("files updated")

#data_array = data.to_numpy(dtype=float)
dropzero = data[np.all(data == 0, axis=1)].index
data2 = data.drop(dropzero, 'index') #drop genes that are not detected in any cells
pvals = pvals.drop(dropzero, 'index') #drop genes that are not detected in any cells, some genes pvals are zero, with non zero counts, why?

coeffs = variation(data2, axis=1)
genes_pcell = np.count_nonzero(data2,0)
avg_oall = np.mean(data2, 1)

data2[data2 == 0] = np.nan
avg_oall[avg_oall == 0] = np.nan
coeffv0s = variation(data2, axis=1, nan_policy='omit')
avg_level = np.nanmean(data2, 1)
count_depth = np.nansum(data2, 1)  # how many times each mRNA is counted across all cells
counts_pcell = np.nansum(data2, 0)  # reads per cell
rank_counts = ss.rankdata(counts_pcell)
num_cells = data2.count(axis = 1) # number of cells each mRNA is detected in
NT = data2.shape[1]  # total ncells there are 17k genes in hydra! 25k cells
n_cells =str(NT)

Pc = 0.05
maxv = (np.nanmax(avg_oall))
minv = (np.nanmin(avg_oall))
ma = np.logspace(np.log10(minv), np.log10(maxv), 5000)
ms = [x / Pc for x in ma]
Yexp = NT * (1 - (np.power((1 - Pc), ms)))
n_starting_avg = np.round(avg_oall * NT / Pc)
n_i_starting = np.round(n_starting_avg / NT)
Pr_zero = (1 - Pc) ** n_i_starting

#BH correction
pcut = 0.01
rankps = ss.rankdata(pvals, method='min')
bhv = pd.DataFrame((rankps/data.shape[0])*pcut)
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
#sig_genes.to_csv(large_root + "/" + label + "significant_genes.csv")

#create frame with only sig gene values
nonsig = bhv.index[bhv['sigp'] == 0]
trimdata = data.drop(nonsig, 'index') #drop genes that are not detected in any cells
trimdata.to_csv(large_root + "/" + label + "trans_sig_data.csv")

num_sig_pvals = np.sum(bhv['sigp'])
f_sig = num_sig_pvals/pvals.shape[0] #of genes that are expressed at all in data
frac_sig = str(f_sig)
num_g = str(pvals.shape[0])
p_cut = str(maxsig)

#select genes for density vis
kn = KneeLocator(ma, Yexp, curve='concave', direction='increasing')
readsmax = kn.knee
g_high = avg_oall[(avg_oall > readsmax)].index #genes with saturated reads, ubiq genes
c_low = num_cells[(num_cells < 50)].index #genes detected in fewer than 50 cells
droplist = g_high.append(c_low)
frames = [pvals,num_cells,avg_oall]
mid_params = pd.concat(frames,axis=1)
mid_params = mid_params.drop(droplist, 'index')
mid_params.columns = ['Pval', 'numcells', 'avg_oall']
mid_data = data2.drop(droplist, 'index')
top10 = mid_params.nsmallest(dnum, 'Pval')
top10_list = list(top10.index.values)
mid_list = list(mid_params.index.values)
min_size = mid_data.shape[0]
if min_size >= 5000:
    step = 0.1
elif min_size >= 500:
    step = .5
elif min_size >= 50:
    step = 10
else:
    step = 50

#gene1 = "CD79B" #far
#gene2 = "SERF2" #middle region
#gene3 = "NDUFA13" #line

#gene1 = "CST3" #far
#gene2 = "CYBA" #middle region
#gene3 = "TMEM59" #line

#gene1 = "FCER1G" #far
#gene2 = "OAZ1" #middle region
#gene3 = "SSU72" #line

gene1_list = ["CD247","FCGR3A","FGFBP2","GNLY","JAK1","HLA-DPA1","CD79B"] #far
gene2_list = ["IL32","S100A6","S100A4","CD52","HCST","GZMM","OAZ1"] #middle region
gene3_list = ["FIS1","JTB","PSMD9","STUB1","TMEM160","FKBP8","SSU72"] #uniform

for gene1, gene2, gene3 in zip(gene1_list,gene2_list,gene3_list):
    # density plots for pairs of deviant and normie genes
    fs = 12  # fontsize
    fig = plt.figure(figsize=(12, 8))
    ax4 = fig.add_subplot(231)
    pv = pvals.loc[gene1, :].values
    pva = str(pv)
    ax4.set_title(gene1 + " p-val: " + pva, fontsize=fs)
    ax4.set_ylabel("freq")
    ax4.set_xlabel("mRNA count across " + label)
    dens1 = cellgroup1.loc[gene1, :]
    dens2 = cellgroup2.loc[gene1, :]
    dens3 = cellgroup3.loc[gene1, :]
    dens4 = cellgroup4.loc[gene1, :]
    dens5 = cellgroup5.loc[gene1, :]
    dens6 = cellgroup6.loc[gene1, :]
    dens7 = cellgroup7.loc[gene1, :]
    dens8 = cellgroup8.loc[gene1, :]
    dens9 = cellgroup9.loc[gene1, :]
    #ax4.set_ylim(0, 0.4)
    #ax4.set_xlim(-2, 100)
    sns.kdeplot(np.array(dens1), bw=0.5, color = 'darkgreen', shade= True)
    sns.kdeplot(np.array(dens2), bw=0.5, color = 'seagreen', shade= True)
    sns.kdeplot(np.array(dens3), bw=0.5, color = 'mediumspringgreen', shade= True)
    sns.kdeplot(np.array(dens4), bw=0.5, color = 'paleturquoise', shade= True)
    sns.kdeplot(np.array(dens5), bw=0.5, color = 'olivedrab', shade= True)
    sns.kdeplot(np.array(dens6), bw=0.5, color = 'lawngreen', shade= True)
    sns.kdeplot(np.array(dens7), bw=0.5, color = 'blue', shade= True)
    sns.kdeplot(np.array(dens8), bw=0.5, color = 'cornflowerblue', shade= True)
    sns.kdeplot(np.array(dens9), bw=0.5, color = 'purple', shade= True)
    plt.legend(labels=['B-cells', 'naive cytotoxic', 'monocytes', 'reg T', 'helper T','NK', 'memory T', 'naive T', 'cytotox T'], loc='upper right')
    plt.tight_layout()

    ax4 = fig.add_subplot(232)
    pv = pvals.loc[gene2, :].values
    pva = str(pv)
    ax4.set_title(gene2 + " p-val: " + pva, fontsize=fs)
    ax4.set_ylabel("freq")
    ax4.set_xlabel("mRNA count across " + label)
    dens1 = cellgroup1.loc[gene2, :]
    dens2 = cellgroup2.loc[gene2, :]
    dens3 = cellgroup3.loc[gene2, :]
    dens4 = cellgroup4.loc[gene2, :]
    dens5 = cellgroup5.loc[gene2, :]
    dens6 = cellgroup6.loc[gene2, :]
    dens7 = cellgroup7.loc[gene2, :]
    dens8 = cellgroup8.loc[gene2, :]
    dens9 = cellgroup9.loc[gene2, :]
    #ax4.set_ylim(0, 0.4)
    #ax4.set_xlim(-2, 100)
    sns.kdeplot(np.array(dens1), bw=0.5, color = 'darkgreen', shade= True)
    sns.kdeplot(np.array(dens2), bw=0.5, color = 'seagreen', shade= True)
    sns.kdeplot(np.array(dens3), bw=0.5, color = 'mediumspringgreen', shade= True)
    sns.kdeplot(np.array(dens4), bw=0.5, color = 'paleturquoise', shade= True)
    sns.kdeplot(np.array(dens5), bw=0.5, color = 'olivedrab', shade= True)
    sns.kdeplot(np.array(dens6), bw=0.5, color = 'lawngreen', shade= True)
    sns.kdeplot(np.array(dens7), bw=0.5, color = 'blue', shade= True)
    sns.kdeplot(np.array(dens8), bw=0.5, color = 'cornflowerblue', shade= True)
    sns.kdeplot(np.array(dens9), bw=0.5, color = 'purple', shade= True)
    plt.legend(labels=['B-cells', 'naive cytotoxic', 'monocytes', 'reg T', 'helper T','NK', 'memory T', 'naive T', 'cytotox T'], loc='upper right')
    plt.tight_layout()

    ax4 = fig.add_subplot(233)
    pv = pvals.loc[gene3, :].values
    pva = str(pv)
    ax4.set_title(gene3 + " p-val: " + pva, fontsize=fs)
    ax4.set_ylabel("freq")
    ax4.set_xlabel("mRNA count across " + label)
    dens1 = cellgroup1.loc[gene3, :]
    dens2 = cellgroup2.loc[gene3, :]
    dens3 = cellgroup3.loc[gene3, :]
    dens4 = cellgroup4.loc[gene3, :]
    dens5 = cellgroup5.loc[gene3, :]
    dens6 = cellgroup6.loc[gene3, :]
    dens7 = cellgroup7.loc[gene3, :]
    dens8 = cellgroup8.loc[gene3, :]
    dens9 = cellgroup9.loc[gene3, :]
    # ax4.set_ylim(0, 0.4)
    # ax4.set_xlim(-2, 100)
    sns.kdeplot(np.array(dens1), bw=0.5, color='darkviolet', shade=True)
    sns.kdeplot(np.array(dens2), bw=0.5, color='m', shade=True)
    sns.kdeplot(np.array(dens3), bw=0.5, color='hotpink', shade=True)
    sns.kdeplot(np.array(dens4), bw=0.5, color='thistle', shade=True)
    sns.kdeplot(np.array(dens5), bw=0.5, color='maroon', shade=True)
    sns.kdeplot(np.array(dens6), bw=0.5, color='firebrick', shade=True)
    sns.kdeplot(np.array(dens7), bw=0.5, color='red', shade=True)
    sns.kdeplot(np.array(dens8), bw=0.5, color='lightcoral', shade=True)
    sns.kdeplot(np.array(dens9), bw=0.5, color='orange', shade=True)
    plt.legend(labels=['B-cells', 'naive cytotoxic', 'monocytes', 'reg T', 'helper T', 'NK', 'memory T', 'naive T',
                       'cytotox T'],
               loc='upper right')
    plt.tight_layout()

    ax3 = fig.add_subplot(234)
    ax3.set_title(label, fontsize=fs)
    ax3.set_ylabel("number of cells")
    ax3.set_xlabel("average mRNA count over all cells")
    ax3.loglog()
    ax3 = plt.plot(ma, Yexp, linewidth=3, color='silver')
    ax3 = plt.scatter(avg_oall, num_cells, color='k', s=100, edgecolors='None', marker=".")
    ax3 = plt.scatter(mid_params.at[gene1,'avg_oall'],mid_params.at[gene1,'numcells'], color = 'cyan', s=100, edgecolors='None', marker=".")
    ax3 = plt.scatter(mid_params.at[gene2,'avg_oall'],mid_params.at[gene2,'numcells'], color = 'cornflowerblue', s=100, edgecolors='None', marker=".")
    ax3 = plt.scatter(mid_params.at[gene3,'avg_oall'],mid_params.at[gene3,'numcells'], color = 'magenta', s=100, edgecolors='None', marker=".")
    plt.tight_layout()
    plt.ylim(1000, 45000)
    plt.xlim(.1, 10)
    plt.savefig(large_root + "/" + label + gene1 + " R01scatter.png")
    plt.savefig(large_root + "/" + label + gene1 + " R01scatter.eps")
