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
sc.settings.verbosity = 0

def rename(A):
    # for gene list A return gene names with spaces instead of - or . characters
    newlist = []
    for gene in A:
        newgene = str(gene)
        a = newgene.replace('.', '-') #R changes - to .
        newlist.append(a)
    all_headers = np.array(newlist)
    return (all_headers)

def PCA_Elbow_fit(data):
    model = PCA().fit(data)
    explained_variance = model.explained_variance_ratio_
    pcs = list(range(1,explained_variance.shape[0]+1))#enumerate(explained_variance,1)
    klm = kl(pcs, explained_variance, S=1.0, curve='convex', direction='decreasing')
    pcs_used = klm.knee
    pc_list = list(range(1,pcs_used+1))
    new_data = PCA(n_components= pcs_used,svd_solver= "arpack").fit_transform(data)
    return pcs_used, new_data, pcs, explained_variance,pc_list


large_root = r"/home/breanne/ddgs/ari"
label = 'testzheng5k_plusmarkers'

#test_cols = [*range(0, 500, 1)]
total_data = pd.read_csv(r"/home/breanne/ddgs/zheng5k_dropzeros.csv", index_col=0) #gene x cells
total_data = total_data.T
markers = pd.read_csv(r"/home/breanne/ddgs/zhengnine/zheng_Panglao_markergenes.csv", index_col=1) #gene list
#cpm_data = pd.read_csv(r"/media/timothyhamilton/My Passport1/Tim_Hamilton/Zheng_FACS_pbmcs/zheng9_5k/zheng5k_cpmlog.csv", index_col=0) #cell x genes? use rows
#WC = pd.read_csv(r"/media/timothyhamilton/My Passport1/Tim_Hamilton/Zheng_FACS_pbmcs/zheng9_5k/zheng9_5k_WC.csv", index_col=0, usecols = [0,1]) #g x c
WC = pd.read_csv(r"/home/breanne/ddgs/zheng5k_oldWilcox_bf8k.csv", index_col=0, nrows=1) #c x g
HVG = pd.read_csv(r"/media/timothyhamilton/My Passport2/Tim_Hamilton/Zheng_FACS_pbmcs/zheng9_5k/zheng5k_HVGs.csv", index_col=0,skiprows=[1156],usecols = [0,1]) #remoe the extra indeces that come up from the duplicate gene
# 'Unnamed: 1156' - is index of 'NA-1' - gene, which is probably the duplicate total data gene that has same name when - are converted to .
DDG = pd.read_csv(r"/media/timothyhamilton/My Passport2/Tim_Hamilton/Zheng_FACS_pbmcs/zheng9_5k/zheng9_5k_DDGs.csv", index_col=0, usecols = [0,1])# g x c
townes = pd.read_csv(r"/media/timothyhamilton/My Passport2/Tim_Hamilton/Zheng_FACS_pbmcs/zheng9_5k/zheng5k_dropzeros_SCRY.csv",index_col = 0,nrows=1)
mdrop = pd.read_csv(r"/media/timothyhamilton/My Passport2/Tim_Hamilton/Zheng_FACS_pbmcs/zheng9_5k/zheng5k_dropzeros_M3Drop.csv",index_col = 0,nrows=1)
nbdrop = pd.read_csv(r"/media/timothyhamilton/My Passport2/Tim_Hamilton/Zheng_FACS_pbmcs/zheng9_5k/zheng5k_dropzeros_NBUMI.csv", index_col=0, nrows = 1)#c x g
nbdisp = pd.read_csv(r"/media/timothyhamilton/My Passport2/Tim_Hamilton/Zheng_FACS_pbmcs/zheng9_5k/zheng5k_dropzeros_NBUMI_High_Var.csv", index_col=0, nrows = 1) # c x g

marker_first = rename(markers.index.tolist())
WC_first= rename(WC.columns.tolist())
all_set= rename(total_data.columns.tolist())
DDG_set= rename(DDG.index.tolist())
M3Drop_set= rename(mdrop.columns.tolist())
Townes_set= rename(townes.columns.tolist())
HVG_first = rename(HVG.index.tolist())
nbdrop_set= rename(nbdrop.columns.tolist())
nbdisp_set = rename(nbdisp.columns.tolist())

HVG_set = [x for x in HVG_first if x in all_set]
marker_set = [x for x in marker_first if x in all_set]
WC_set = [x for x in WC_first if x in all_set]
#print(all_set)
#print(marker_set)
marker_set = list(marker_set)
WC_set= list(WC_set) + marker_set
DDG_set= list(DDG_set) + marker_set
M3Drop_set= list(M3Drop_set) + marker_set
Townes_set= list(Townes_set) + marker_set
HVG_set = list(HVG_set) + marker_set
nbdrop_set= list(nbdrop_set) + marker_set
nbdisp_set = list(nbdisp_set) + marker_set

title_list = ["marker genes","WC","DDGs","HVGs","M3Drop","Townes","NBdrop","NBdisp"]
gene_list = [marker_set, WC_set, DDG_set, HVG_set, M3Drop_set, Townes_set, nbdrop_set, nbdisp_set]

for gene_group, glabel in zip(gene_list, title_list):
    print(glabel)
    print(len(gene_group))

barcodes = np.array(total_data.index)
genenames = all_set
data = total_data.to_numpy(dtype=int)
data = pd.DataFrame(data=data, index=barcodes, columns=genenames)
adata = sc.AnnData(data)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

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
total_data= pd.DataFrame(total_data, index = facs_headers, columns = genenames)
cpm_data = pd.DataFrame(data=adata.X, index = facs_headers, columns=genenames)
print(cpm_data.head())

step = 1
reslist = list(range(1, 11, step))
LV_lables = [str(x) for x in reslist]
LV_lables2 = ['p'+x for x in LV_lables]

data_list = [total_data, cpm_data]
countname = ['_raw_counts_', '_cpm_log_counts_']

print(total_data.shape)
print(cpm_data.shape)

#cycle through different normalizations
for data_group, cname in zip(data_list,countname):

    # cycle through different feature sets
    for glabel, gname in zip(gene_list,title_list):
        test_group = pd.Series(glabel)
        new_data = data_group.copy()
        new_data = new_data.loc[:, test_group]

        #PCA transform
        new_data = new_data.fillna(0) #there shouldnt be any though
        nm = new_data.values
        dim, new_matrix, pc_ax, pc_ay, col_labels = PCA_Elbow_fit(nm)
        columns = ["PC_" + str(i) for i in col_labels]
        output_path = large_root + gname + cname + "_PCA_" + str(dim) + ".csv"
        new_frame = pd.DataFrame(new_matrix, index=new_data.index.values.tolist(), columns=columns)
        print(new_frame.head())
        print(new_frame.shape)
        # data=data.transpose()
        new_frame.to_csv(output_path)
        del new_frame
        del nm

        #do louvain clustering for each res in reslist on high D and PCA reduced data at elbow dim
        for res,reslabel,reslabel2 in zip(reslist,LV_lables,LV_lables2):
            res = res/10
            print("Calculating Clustering")
            keys_list = []
            raw_clustering = []
            PC_clustering = []
            data = sc.AnnData(new_data)
            sc.pp.neighbors(data, n_neighbors=20, use_rep='X')  # this X is louvain on non PCA data
            sc.tl.louvain(data, resolution=res, key_added=reslabel)
            keys_list.append(reslabel)
            raw_clustering = data.obs[reslabel]
            rawc = raw_clustering.tolist()
            top_clusterin[reslabel] = rawc

            sc.pp.neighbors(data, n_neighbors=20, n_pcs=dim)  # this X is louvain on PCA transform
            sc.tl.louvain(data, resolution=res, key_added=reslabel2)
            keys_list.append(reslabel2)
            PC_clustering = data.obs[reslabel2]
            PCc = PC_clustering.tolist()
            PC_clusterin[reslabel2] = PCc

            data.obs['facs'] = top_clusterin.FACS
            sc.tl.pca(data, n_comps=50)
            sc.tl.tsne(data, n_pcs=50)

            # make eps files I want with colors by avoiding scanpy
            tsne_cords = pd.DataFrame(data.obsm["X_tsne"], index=facs_headers)
            tsne_cords.rename(columns={tsne_cords.columns[0]: "tsne1"}, inplace=True)
            tsne_cords.rename(columns={tsne_cords.columns[1]: "tsne2"}, inplace=True)
            tsne_cords['LV'] = top_clusterin[reslabel]
            tsne_cords['LV-PC'] = PC_clusterin[reslabel2]
            tsne_cords['FACS'] = top_clusterin['FACS']

            # there are 9 LV, 5 FACS ############ will have to edit this
            LV_colors = ['dodgerblue', 'darkgreen', 'firebrick', 'rebeccapurple', \
                           'cyan', 'grey', 'darkorange', 'magenta']
            FACS_colors = ['dodgerblue', 'darkgreen', 'firebrick', 'rebeccapurple', \
                           'cyan', 'grey', 'darkorange', 'magenta']
            fs = 10
            sns.set_palette(sns.color_palette(LV_colors))
            sns.lmplot(x="tsne1", y="tsne2", data=tsne_cords, fit_reg=False, hue='LV', legend=None, scatter_kws={"s": 2})
            ax = plt.gca()
            ax.set_title(label + gname + cname + "\n" + "colored by LV res 0." + reslabel, fontsize=fs)
            ax.set_ylabel("tSNE_2")
            ax.set_xlabel("tSNE_1")
            plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
            plt.tight_layout()
            sres = res*10
            #plt.savefig(large_root + "/" + label + gname + cname + reslabel + "LV_TSNE.png")
            #plt.savefig(large_root + "/" + label + gname + cname + reslabel + "LV_TSNE.eps")

            sns.set_palette(sns.color_palette(LV_colors))
            sns.lmplot(x="tsne1", y="tsne2", data=tsne_cords, fit_reg=False, hue='LV-PC', legend=None, scatter_kws={"s": 2})
            ax = plt.gca()
            ax.set_title(label + "PCA-" + gname + cname + "\n" + " colored by LV res 0." + reslabel, fontsize=fs)
            ax.set_ylabel("tSNE_2")
            ax.set_xlabel("tSNE_1")
            plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
            plt.tight_layout()
            #plt.savefig(large_root + "/" + label + "PCA-" + gname + cname + reslabel + "LV_TSNE.png")
            #plt.savefig(large_root + "/" + label + "PCA-" + gname + cname + reslabel + "LV_TSNE.eps")

        sns.set_palette(sns.color_palette(FACS_colors))
        sns.lmplot(x="tsne1", y="tsne2", data=tsne_cords, fit_reg=False, hue='FACS', legend=None, scatter_kws={"s": 2})
        ax = plt.gca()
        ax.set_title(label + gname + cname + " colored byFACS ", fontsize=fs)
        ax.set_ylabel("tSNE_2")
        ax.set_xlabel("tSNE_1")
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.tight_layout()
        #plt.savefig(large_root + "/" + label + gname + cname + "\n" + "FACS_TSNE.png")
        #plt.savefig(large_root + "/" + label + gname + cname + "\n" + "FACS_TSNE.eps")

        xxlabels = list(top_clusterin.columns)
        xlabels = xxlabels[1:]
        #PC_clustering = data.obs[LVres]
        adjusted_rand = []
        adjusted_randPC = []
        for reso, reslabel, reslabel2 in zip(reslist, LV_lables, LV_lables2):
            trial = adjusted_rand_score(top_clusterin[reslabel], top_clusterin.FACS)
            PCtrial = adjusted_rand_score(PC_clusterin[reslabel2], top_clusterin.FACS)
            adjusted_randPC.append(PCtrial)
            adjusted_rand.append(trial)
        randdata = pd.DataFrame({'LVres': xlabels, 'ARI': adjusted_rand})
        randdata.to_csv(large_root + "/" + label + gname + cname + "ARIs.csv")
        randdataPC = pd.DataFrame({'LVres': xlabels, 'ARI': adjusted_randPC})
        randdataPC.to_csv(large_root + "/" + label + gname + cname + "ARIsPC.csv")

        fs = 12
        fig = plt.figure(figsize=(5, 8))
        ax = fig.add_subplot(111)
        ax = randdata.plot.bar(x='LVres', y='ARI', rot=0, width=0.75)
        ax.set_title("ARI: FACs labels vs Clustering Resolution \n" + label + gname + cname, fontsize=fs)
        ax.set_ylabel("Adjusted Rand Index")
        ax.set_xlabel("LV resolution")
        ax.set_ylim(0, 1)
        plt.grid()
        #plt.savefig(large_root + "/" + label + gname + cname + "resARI.png")
        plt.savefig(large_root + "/" + label + gname + cname + "resARI.eps")

        fs = 12
        fig = plt.figure(figsize=(5, 8))
        ax = fig.add_subplot(111)
        ax = randdataPC.plot.bar(x='LVres', y='ARI', rot=0, width=0.75)
        ax.set_title("ARI: FACs labels vs Clustering Resolution \n" + "PC data" + label + gname + cname, fontsize=fs)
        ax.set_ylabel("Adjusted Rand Index")
        ax.set_xlabel("LV resolution")
        ax.set_ylim(0, 1)
        plt.grid()
        #plt.savefig(large_root + "/" + label + gname + cname + "PCresARI.png")
        plt.savefig(large_root + "/" + label + gname + cname + "PCresARI.eps")