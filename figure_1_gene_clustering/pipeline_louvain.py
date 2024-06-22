import pandas as pd
import numpy as np
import scanpy as sc
from sklearn.metrics.cluster import adjusted_rand_score
import matplotlib.pyplot as plt
import sys
import os
import seaborn as sns

# read the data from 10x .mtx:
def louvain_subclustering(sourcepath,source,delim):
    large_root=source
    print ("Doing: "+ source)
    newpath1 = sourcepath+"/subclusters/"
    newpath2 = sourcepath+"/figures/"
    print("making new directories")
    try:
        print("making directory for random subcluster")
        os.makedirs(newpath1)
        os.makedirs(newpath2)
    except OSError:
        print("Directories for subclusters and figures already exist, moving on ...")
    sc.settings.autoshow=False
    sc.settings.figdir = newpath2
    print("Reading data...")
        #marker_genes = ["MYCL","EPCAM","REST", "NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4","HES1","CD44","ASCL1","DDC","GRP","NEUROD1","NEUROD2","NEUROD4", "MYC","CDH5","CLDN4","POU2F3","EZH2"]
    adata= sc.read_csv(source, first_column_names = True, delimiter = delim)
    adata = adata.T
    print("Calculating Neighbors")
    sc.pp.neighbors(adata,n_neighbors = 20,use_rep = 'X')
    sc.tl.louvain(adata,resolution = 0.5)
        #new_data = new_data[new_data.obs['louvain'].isin(['0','1','2', '3'])]
        #print(new_data.layers)
    print("Removing Subclusters")
    qol_indices=[] ##6/19 newline
    qol_cols=['Cell Count', 'Gene Count'] ##6/19 newline
    qol_list=[] ##6/19 newline
    clust_list = sorted([int(i) for i in list(set(adata.obs["louvain"]))])
    for id in clust_list:
        sub_data = adata[adata.obs['louvain'].isin([str(id)])]
        qol_list.append([sub_data.X.shape[0],sub_data.X.shape[1]]) ##6/19 newline
        qol_indices.append(str(id)) ##6/19 newline
        sub_frame = pd.DataFrame(data=sub_data.X, index=sub_data.obs_names, columns=sub_data.var_names)
        sub_frame.to_csv(newpath1+"subcluster_"+str(id)+".csv")
    qol_np=np.asarray(qol_list) ##6/19 newline
    qol_df=pd.DataFrame(qol_np,index=qol_indices,columns=qol_cols) ##6/19 newline
    qol_df.to_csv(newpath2+'subcluster_data.csv') ##6/19 newline
    fig,ax=plt.subplots(figsize=(20,4)) ##6/19 newline
    fig.patch.set_visible(False) ##6/19 newline
    ax.xaxis.set_visible(False) ##6/19 newline
    ax.yaxis.set_visible(False) ##6/19 newline
    ax.table(cellText=qol_df.values,colLabels=qol_df.columns,rowLabels=qol_df.index,loc='center') ##6/19 newline
    plt.title("Summary of Samples") ##6/19 newline
    fig.savefig(newpath2+'subcluster_data.png',dpi=200) ##6/19 newline
    plt.close() ##6/19 newline

    adata.obs['avg_exp'] = adata.X.mean(1)  # avg counts per gene
    adata.obs['log_avg_exp'] = np.log(adata.X.mean(1))  # avg counts per gene
    #t1 = sc.pl.violin(adata, 'avg_exp', groupby='louvain', size=2, log=True, cut=0,save='avg_count.png')  ##6/19 newline #
    #t1 = sc.pl.violin(adata, 'avg_exp', groupby='louvain', size=2, log=True, cut=0,save='avg_count.eps')  ##6/19 newline #
    t1 = sc.pl.violin(adata, 'log_avg_exp', groupby='louvain', size=2, log=True, save='log_avg_count.png')  ##6/19 newline #
    t1 = sc.pl.violin(adata, 'log_avg_exp', groupby='louvain', size=2, log=True, save='log_avg_count.eps')  ##6/19 newline #

    adata.obs['n_counts'] = adata.X.sum(1)  ##6/19 newline # how many counts each gene has
    adata.obs['log_counts'] = np.log(adata.obs['n_counts'])  ##6/19 newline
    adata.obs['n_genes'] = (adata.X > 0).sum(1)  ##6/19 newline # num of cells each gene is exp in
    # n cells
    adata.obs['frac_cells'] = (adata.X > 0).sum(1) / adata.X.shape[1]
    adata.obs['log_frac_cells'] = np.log((adata.X > 0).sum(1) / adata.X.shape[1])
    t1 = sc.pl.violin(adata, 'log_frac_cells', groupby='louvain', size=2, log=True, save='log_frac_cells.png')  ##6/19 newline #
    t1 = sc.pl.violin(adata, 'log_frac_cells', groupby='louvain', size=2, log=True, save='log_frac_cells.eps')  ##6/19 newline #
    #t1 = sc.pl.violin(adata, 'frac_cells', groupby='louvain', size=2, log=True, cut=0, save='frac_cells.png')  ##6/19 newline #
    #t1 = sc.pl.violin(adata, 'frac_cells', groupby='louvain', size=2, log=True, cut=0, save='frac_cells.eps')  ##6/19 newline #

    mt_gene_mask = [gene.startswith('mt-') for gene in adata.var_names] ##6/19 newline
    adata.obs['mt_frac'] = adata.X[:,mt_gene_mask].sum(1) / adata.obs['n_counts'] ##6/19 newline
    t1=sc.pl.violin(adata,'n_counts',groupby='louvain',size=2,log=True,cut=0,save='counts_p_gene.png') ##6/19 newline #
    t1=sc.pl.violin(adata,'n_counts',groupby='louvain',size=2,log=True,cut=0,save='counts_p_gene.eps') ##6/19 newline #

    t2=sc.pl.violin(adata,'mt_frac',groupby='louvain',save='q2plot.png') ##6/19 newline
    p1=sc.pl.scatter(adata,'n_counts','n_genes',color='mt_frac',save='q3plot.png') ##6/19 newline
    p2=sc.pl.scatter(adata[adata.obs['n_counts']<10000],'n_counts','n_genes',color='mt_frac',save='q4plot.png') ##6/19 newline
    p3=sns.distplot(adata.obs['n_counts'],kde=False).get_figure() ##6/19 newline
    p3.savefig(newpath2+'counts1.png') ##6/19 newline
    p3.clf() ##6/19 newline
    p4=sns.distplot(adata.obs['n_counts'][adata.obs['n_counts']<4000],kde=False,bins=60).get_figure() ##6/19 newline
    p4.savefig(newpath2+'counts2.png') ##6/19 newline
    p4.clf() ##6/19 newline
    p5=sns.distplot(adata.obs['n_counts'][adata.obs['n_counts']>10000],kde=False,bins=60).get_figure() ##6/19 newline
    p5.savefig(newpath2+'counts3.png') ##6/19 newline
    p5.clf() ##6/19 newline
    p6=sns.distplot(adata.obs['n_genes'],kde=False,bins=60).get_figure() ##6/19 newline
    p6.savefig(newpath2+'genecounts1.png') ##6/19 newline
    p6.clf() ##6/19 newline
    p7=sns.distplot(adata.obs['n_genes'][adata.obs['n_genes']<1000],kde=False,bins=60).get_figure() ##6/19 newline
    p7.savefig(newpath2+'genecounts2.png') ##6/19 newline
    p7.clf() ##6/19 newline


    sc.pp.normalize_total(adata,target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata
    sc.pp.highly_variable_genes(adata,subset = True, n_top_genes = 2000)
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps= 50)
    sc.tl.tsne(adata,n_pcs=50)
    print("making tsne")
    sc.pl.tsne(adata, color = ['louvain'], save= "louvain_tsne_plot.png")
    print("making umap") ##6/19 newline
    sc.tl.umap(adata) ##6/19 newline
    sc.pl.umap(adata, color=['louvain'], save="louvain_umap_plot.png") ##6/19 newline
    print("making paga")
    sc.tl.paga(adata) ##6/19 newline
    sc.pl.paga(adata,color=['louvain'],save="louvain_paga_plot.png") ##6/19 newline
    print("making heatmap")
    #sc.pl.heatmap(adata,marker_genes,groupby='louvain',save = "_"+".png",standard_scale='var')
    print("Making violin plot")
   # sc.pl.stacked_violin(adata,marker_genes,groupby='louvain',save = "_"+".png",standard_scale='var')
    sc.tl.rank_genes_groups(adata, groupby='louvain')
    sc.tl.filter_rank_genes_groups(adata)
    sc.pl.rank_genes_groups(adata,key='rank_genes_groups_filtered',save = "_"+"_marker_gene.png")
    sc.pl.rank_genes_groups_heatmap(adata,key='rank_genes_groups_filtered',save = "_"+"_marker_gene.png",groupby='louvain',standard_scale='var')
    sc.pl.rank_genes_groups_stacked_violin(adata,key='rank_genes_groups_filtered',save = "_"+"_marker_gene.png",groupby='louvain',standard_scale='var')

    sc.pl.rank_genes_groups_heatmap(adata,key='rank_genes_groups_filtered',save = "_"+"_marker_gene.eps",groupby='louvain',standard_scale='var',cmap = "rocket")
    sc.pl.rank_genes_groups_stacked_violin(adata,key='rank_genes_groups_filtered',save = "_"+"_marker_gene.eps",groupby='louvain',standard_scale='var', cmap = "rocket")


### if you want to read a loom file:
# adata = sc.read_loom(filename)

# Do a nearest neighbor search:
