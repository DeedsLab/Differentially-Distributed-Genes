import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.colors
import scipy.stats as ss
import seaborn as sns
from kneed import KneeLocator
from scipy.stats import variation

large_root = r"C:\Users\bsparta\Desktop\data\clustering_validation"
label = 'Whole_hydra'

cell_types = pd.read_csv(r"C:\Users\bsparta\Desktop\data\R01\rawdata\Hydra_Raw_Cell_Types.csv", usecols=[0,1,2,3] )
cell_type_array = cell_types['Cell_Type']
total_data = pd.read_csv(r"C:\Users\bsparta\Desktop\data\R01\rawdata\Hydra_full_RawTrans.csv")
sig_data = pd.read_csv(r"C:\Users\bsparta\Desktop\data\R01\sig_gene_lists_0_01_bh_total\Hydratrans_sig_data.csv")
HVG_list = pd.read_csv(r"C:\Users\bsparta\Desktop\data\R01\Seurat_results\Hydra_Var_features.csv")
pval = pd.read_csv(r"C:\Users\bsparta\Desktop\data\R01\pval_fast\hydra_whole24458_fast_Pvalues.csv", header = None, names = ["gene", "Pval"])

total_data = total_data.set_index('Unnamed: 0')
cell_headers = np.array(total_data.columns)

all_genes = total_data.index.tolist()
newlist = []
for gene in all_genes:
    newgene = str(gene)
    a = newgene.replace('.','').replace('-','')
    newlist.append(a)
all_headers = np.array(newlist)
Total_data = total_data.values
Total_Data=pd.DataFrame(Total_data, index=all_headers, columns = cell_type_array)

sig_data = sig_data.set_index('Unnamed: 0')
sig_head = sig_data.index.tolist()
newlist = []
for gene in sig_head:
    newgene = str(gene)
    a = newgene.replace('.','').replace('-','')
    newlist.append(a)
sig_headers = np.array(newlist)
ddg_data = sig_data.values
DDG_data=pd.DataFrame(ddg_data, index=sig_headers, columns = cell_type_array)

HVG_list = HVG_list.set_index('Unnamed: 0')
HVG_list=HVG_list.T
hvg_head = HVG_list.index.tolist()
newlist = []
for gene in hvg_head:
    newgene = str(gene)
    a = newgene.replace('.','').replace('-','')
    newlist.append(a)
hvg_headers = np.array(newlist)
hvg_data = HVG_list.values
HVG_data=pd.DataFrame(hvg_data, index=hvg_headers, columns = cell_type_array)

DDG_data = DDG_data.T
HVG_data = HVG_data.T

hvg_headers.tolist()
sig_headers.tolist()
comb_feat = []
for gene in hvg_headers:
    if gene not in sig_headers:
        newgene = gene
        comb_feat.append(newgene)
comb_feat.extend(sig_headers)

print(np.shape(sig_headers)) #ddgs
print(np.shape(hvg_headers)) #hvgs
print(np.shape(comb_feat)) #comb

cell_type_list = []
for cell in cell_type_array:
    if cell not in cell_type_list:
        newcell = cell
        cell_type_list.append(newcell)

total_data = Total_Data
data_list = []
for celltype in cell_type_list:
    celltypename = celltype
    celltype = total_data[celltypename]
    data_list.append(celltype)

data_list.append(total_data)
cell_type_list.append('all_cells')

pval = pval.drop(columns=['gene'])
pval.index = all_headers

print("files updated")

for cell_group, clabel in zip(data_list,cell_type_list):

    ddg_data = cell_group.filter(sig_headers, axis=0)
    hvg_data = cell_group.filter(hvg_headers, axis=0)
    all_data = cell_group.filter(all_headers, axis=0)
    gene_list = [ddg_data, hvg_data, all_data]
    gene_names = ['_DDGs', '_HVGs', '_all_genes']

    for gene_group, glabel in zip(gene_list, gene_names):
        data = gene_group
        pvals = pval
        p_headers = np.array(data.index)
        pvals = pvals.filter(p_headers, axis = 0)

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

        #create frame with only sig gene values
        nonsig = bhv.index[bhv['sigp'] == 0]
        trimdata = data.drop(nonsig, 'index') #drop genes that are not detected in any cells

        num_sig_pvals = np.sum(bhv['sigp'])
        f_sig = num_sig_pvals/pvals.shape[0] #of genes that are expressed at all in data
        frac_sig = str(f_sig)
        num_g = str(pvals.shape[0])
        p_cut = str(maxsig)

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

        #plot summary plots
        cm = plt.cm.get_cmap('cool')
        cm_reverse = cm.reversed()
        fs = 14
        fig = plt.figure(figsize=(24, 8))

        ax = fig.add_subplot(141)
        ax.set_title(label + clabel + glabel, fontsize=fs)
        ax.set_ylabel("counts per cell")
        ax.set_xlabel("rank")
        sc = plt.scatter(rank_counts, counts_pcell, s=30, marker="o")
        plt.tight_layout()

        ax = fig.add_subplot(142)
        ax.set_title(label + clabel + glabel, fontsize=fs)
        ax.set_ylabel("counts per cell")
        ax.set_xlabel("genes per cell")
        sc = plt.scatter(genes_pcell, counts_pcell, s=30, marker="o")
        plt.tight_layout()

        ax = fig.add_subplot(143)
        ax.set_title(label + clabel + glabel + "\n BH corr p-value threshold: " + p_cut + "\nfraction significant: " + frac_sig, fontsize=fs)
        ax.set_ylabel("number of cells")
        ax.set_xlabel("average mRNA count over all cells")
        ax.loglog()
        sc1 = plt.scatter(avg_oall, num_cells, c=pmap2, s=30, edgecolors='',
                          marker=".", cmap=cm, norm=matplotlib.colors.LogNorm())
        plt.colorbar(sc1, label="p-values", orientation='horizontal', extend = 'both')
        sc1 = plt.plot(ma, Yexp, linewidth=3)
        plt.tight_layout()

        ax = fig.add_subplot(144)
        ax.set_title(label + clabel + glabel + "\nnumber genes: " + num_g + "\nnumber cells: " + n_cells, fontsize=fs)
        ax.set_ylabel("number of cells")
        ax.set_xlabel("average mRNA count over all cells")
        ax.loglog()
        sc4 = plt.scatter(avg_oall, num_cells, c=coeffs, vmin=np.nanmin(coeffs), vmax=np.nanmax(coeffs), s=30, edgecolors='',
                          marker=".", cmap=cm_reverse, norm=matplotlib.colors.LogNorm())
        plt.colorbar(sc4, label="coeff variation over all", orientation='horizontal', extend = 'both')
        sc4 = plt.plot(ma, Yexp, linewidth=3)
        plt.tight_layout()

        plt.savefig(large_root + "/" + label + clabel + glabel + "Log countd pval scatter.png")
        plt.savefig(large_root + "/" + label + clabel + glabel + "Log countd pval scatter.eps")
        fig.clear()
