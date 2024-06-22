import scvelo as scv
import pandas as pd
import scanpy as sc
import numpy as np

large_root = r"/home/breannesparta/ddgs/pancreas"
scv.set_figure_params('scvelo')
scv.settings.figdir=large_root

#load example dataset -- maybe use one of theirs - pancreas - calc ddgs / hvgs for these idk
adata = scv.datasets.pancreas()
#adata.write_csvs(large_root, sep=',')
#spliced counts in adata.X
#large_root = r"/Users/breannesparta/Desktop/results/ddgs/pancreas"
#barcodes = pd.read_csv(large_root + "/obs.csv", index_col=0) #gene list
#genes = pd.read_csv(large_root + "/var.csv", index_col=0) #gene list
#bar_list = barcodes.index.tolist()
#gene_list = genes.index.tolist()
#b = adata.layers['spliced']
#pancreas_data = pd.DataFrame(data=b.toarray(), index=bar_list, columns=gene_list)
#pancreas_data.to_csv(large_root + "/pancreas.csv")

#filter and normalize and select genes in this step -- here add ddgs
#scv.pp.filter_and_normalize(adata, **params)


ddg_list = pd.read_csv(large_root + r"/pancreas_DDGs.csv", index_col=0, usecols=[0,1]) #gene x cells
ddg_list = ddg_list.index.tolist()
#ddg_raw = adata[:,ddg_list] this works but we need to normalize and this doenst work after filtering

scv.pp.filter_genes(adata, min_shared_counts=20)
scv.pp.normalize_per_cell(adata)

obs_list = adata.var.index.tolist()
ddg_filt = [x for x in ddg_list if x in obs_list]

bdata = adata[:,ddg_filt]
#scv.pp.filter_genes_dispersion(adata, n_top_genes=2000) # skip this for DDG; filter based on DDGs
scv.pp.log1p(bdata)

#calc mean and variance
scv.pp.moments(bdata, n_pcs=30, n_neighbors=30)

#solve velocity model
#scv.tl.velocity(adata, mode='stochastic', **params)

scv.tl.velocity(bdata, mode='dynamical')
scv.tl.velocity_graph(bdata)

"""The velocities are projected into a lower-dimensional embedding by translating them into likely cell transitions.
 That is, for each velocity vector we find the likely cell transitions that are in accordance with that direction.
  The probabilities of one cell transitioning into another cell are computed using cosine correlation 
  (between the potential cell transition and the velocity vector) and are stored in a matrix denoted as velocity graph:"""

adata.write(large_root + 'hvg_model_pancreas.h5ad', compression='gzip')
#adata = scv.read('data/pancreas.h5ad')
scv.pl.velocity_embedding_stream(bdata, basis='umap', save='ddg_umap_vel.png')

scv.tl.recover_dynamics(bdata)
scv.tl.latent_time(bdata)
scv.pl.scatter(bdata, color='latent_time', color_map='gnuplot', size=80, save='ddg_latenttime.png')

top_genes = bdata.var['fit_likelihood'].sort_values(ascending=False).index[:300]
print(top_genes)
topddg = top_genes.tolist()
topddg_frame = pd.DataFrame(topddg)
topddg_frame.to_csv(large_root + "/pancreas_ddg_top300latenttime.csv")
scv.pl.heatmap(bdata, var_names=top_genes, sortby='latent_time', col_color='clusters', n_convolve=100, save='ddg_topgenes.png')