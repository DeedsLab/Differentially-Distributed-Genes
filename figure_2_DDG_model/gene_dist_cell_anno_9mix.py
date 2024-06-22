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
dnum = 400 #number of density graphs to print

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
data2 = data2.fillna(0) #####################3added this

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
bhv = pd.DataFrame((rankps/data2.shape[0])*pcut) ################## this was bug, it was data before zero genes dropped
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
sig_genes.to_csv(large_root + "/" + label + "DDGs.csv")

#create frame with only sig gene values
nonsig = bhv.index[bhv['sigp'] == 0]
trimdata = data2.drop(nonsig, 'index') ###### also bug used data not data2
trimdata.to_csv(large_root + "/" + label + "DDG_data.csv")

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

for gene in top10_list:
    yval = mid_params.at[gene, 'avg_oall']
    subgroup = mid_params.loc[(mid_params['avg_oall'] >= (yval-step)) & (mid_params['avg_oall'] <= (yval+step))]
    max_ex = subgroup.nlargest(1, 'Pval') # out of list of index find the one with max pval
    ex_list = list(max_ex.index.values)
    # density plots for pairs of deviant and normie genes
    fs = 12  # fontsize
    fig = plt.figure(figsize=(10, 8))
    ax4 = fig.add_subplot(231)
    pv = pvals.loc[gene, :].values
    pva = str(pv)
    ax4.set_title(gene + " p-val: " + pva, fontsize=fs)
    ax4.set_ylabel("freq")
    ax4.set_xlabel("mRNA count across " + label)
    dens = data2.loc[gene, :]
    #ax4.set_ylim(0, 0.4)
    #ax4.set_xlim(-2, 100)
    sns.kdeplot(np.array(dens), bw=0.5, color = 'cyan', shade= True)
    plt.tight_layout()

    ax4 = fig.add_subplot(232)
    pv = pvals.loc[gene, :].values
    pva = str(pv)
    ax4.set_title(gene + " p-val: " + pva, fontsize=fs)
    ax4.set_ylabel("freq")
    ax4.set_xlabel("mRNA count across " + label)
    dens1 = cellgroup1.loc[gene, :]
    dens2 = cellgroup2.loc[gene, :]
    dens3 = cellgroup3.loc[gene, :]
    dens4 = cellgroup4.loc[gene, :]
    dens5 = cellgroup5.loc[gene, :]
    dens6 = cellgroup6.loc[gene, :]
    dens7 = cellgroup7.loc[gene, :]
    dens8 = cellgroup8.loc[gene, :]
    dens9 = cellgroup9.loc[gene, :]
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

    for normie in ex_list:
        ax2 = fig.add_subplot(234)
        pv = pvals.loc[normie, :].values
        pva = str(pv)
        ax2.set_title(normie + " p-val: " + pva, fontsize=fs)
        ax2.set_ylabel("freq")
        ax2.set_xlabel("mRNA count across " + label)
        dens = data2.loc[normie, :]
        #ax2.set_ylim(0, 0.40)
        #ax2.set_xlim(-2, 40)
        sns.kdeplot(np.array(dens), bw=0.5, color = 'magenta', shade = True)
        plt.tight_layout()

        ax4 = fig.add_subplot(235)
        pv = pvals.loc[normie, :].values
        pva = str(pv)
        ax4.set_title(normie + " p-val: " + pva, fontsize=fs)
        ax4.set_ylabel("freq")
        ax4.set_xlabel("mRNA count across " + label)
        dens1 = cellgroup1.loc[normie, :]
        dens2 = cellgroup2.loc[normie, :]
        dens3 = cellgroup3.loc[normie, :]
        dens4 = cellgroup4.loc[normie, :]
        dens5 = cellgroup5.loc[normie, :]
        dens6 = cellgroup6.loc[normie, :]
        dens7 = cellgroup7.loc[normie, :]
        dens8 = cellgroup8.loc[normie, :]
        dens9 = cellgroup9.loc[normie, :]
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
        plt.legend(labels=['B-cells', 'naive cytotoxic', 'monocytes', 'reg T', 'helper T', 'NK', 'memory T', 'naive T', 'cytotox T'],
                   loc='upper right')
        plt.tight_layout()

        ax3 = fig.add_subplot(233)
        ax3.set_title(label, fontsize=fs)
        ax3.set_ylabel("number of cells")
        ax3.set_xlabel("average mRNA count over all cells")
        ax3.loglog()
        sc2 = plt.plot(ma, Yexp, linewidth=3, color='silver')
        sc2 = plt.scatter(avg_oall, num_cells, color='k', s=100, edgecolors='None', marker=".")
        sc2 = plt.scatter(mid_params.at[gene,'avg_oall'],mid_params.at[gene,'numcells'], color = 'cyan', s=100, edgecolors='None', marker=".")
        sc2 = plt.scatter(mid_params.at[normie,'avg_oall'],mid_params.at[normie,'numcells'], color = 'magenta', s=100, edgecolors='None', marker=".")
        plt.tight_layout()
        ax3.set_ylim(50)
        ax3.set_xlim(.1, 50)
    plt.savefig(large_root + "/" + label + gene + " density.png")
    plt.savefig(large_root + "/" + label + gene + " density.eps")
    fig.clear()

#plot scatters of data
#replace pvalues of zero with very small number
pp = pvals.to_numpy(dtype=float)
pmap = pp.tolist()
pnozero = pvals[~np.all(pvals == 0, axis=1)] #bec some pvals are still zero..
zz = pnozero.to_numpy(dtype=float)
z = zz.tolist()
minz = list(z)
minz.remove(np.nanmin(minz))
while np.min(minz) < 1.79769313486e-308: #float limit
    minz.remove(np.nanmin(minz))
ppp = pvals.to_numpy(dtype=object)
pmap2 = np.concatenate(ppp)
pmap2[pmap2 < np.nanmin(minz)] = np.nanmin(minz)

# plot summary plots
cm = plt.cm.get_cmap('cool')
cm_reverse = cm.reversed()
fs = 14
fig = plt.figure(figsize=(24, 8))

ax = fig.add_subplot(141)
ax.set_title(label, fontsize=fs)
ax.set_ylabel("counts per cell")
ax.set_xlabel("rank")
sc = plt.scatter(rank_counts, counts_pcell, s=30, marker="o")
plt.tight_layout()

ax = fig.add_subplot(142)
ax.set_title(label, fontsize=fs)
ax.set_ylabel("counts per cell")
ax.set_xlabel("genes per cell")
sc = plt.scatter(genes_pcell, counts_pcell, s=30, marker="o")
plt.tight_layout()

ax = fig.add_subplot(143)
ax.set_title(label + "\n BH corr p-value threshold: " + p_cut + "\nfraction significant: " + frac_sig,
             fontsize=fs)
ax.set_ylabel("number of cells")
ax.set_xlabel("average mRNA count over all cells")
ax.loglog()
sc1 = plt.scatter(avg_oall, num_cells, c=pmap2, s=70, edgecolors='None',
                  marker=".", cmap=cm, norm=matplotlib.colors.LogNorm())
#sc1 = plt.autoscale(enable=True, axis='both', tight=True)
plt.colorbar(sc1, label="p-values", orientation='horizontal', extend='both')
sc1 = plt.plot(ma, Yexp, linewidth=3)
plt.tight_layout()

ax = fig.add_subplot(144)
ax.set_title(label + "\nnumber genes: " + num_g + "\nnumber cells: " + n_cells, fontsize=fs)
ax.set_ylabel("number of cells")
ax.set_xlabel("average mRNA count over all cells")
ax.loglog()
sc4 = plt.scatter(avg_oall, num_cells, c=coeffs, vmin=np.nanmin(coeffs), vmax=np.nanmax(coeffs), s=70, edgecolors='None',
                  marker=".", cmap=cm_reverse, norm=matplotlib.colors.LogNorm())
#sc4 = plt.autoscale(enable=True, axis='both', tight=True)
plt.colorbar(sc4, label="coeff variation over all", orientation='horizontal', extend='both')
sc4 = plt.plot(ma, Yexp, linewidth=3)
plt.tight_layout()

plt.savefig(large_root + "/" + label + "Log countd pval scatter.png")
plt.savefig(large_root + "/" + label + "Log countd pval scatter.eps")
fig.clear()