import scvelo as scv
import pandas as pd
import scanpy as sc

large_root = r"/home/breannesparta/ddgs/dentategyrus"
scv.set_figure_params('scvelo')
scv.settings.figdir=large_root

#load example dataset -- maybe use one of theirs - pancreas - calc ddgs / hvgs for these idk
adata = scv.datasets.dentategyrus()
adata.write_csvs(large_root, sep=',')
#spliced counts in adata.X
#large_root = r"/Users/breannesparta/Desktop/results/ddgs/pancreas"
barcodes = pd.read_csv(large_root + "/obs.csv", index_col=0) #gene list
genes = pd.read_csv(large_root + "/var.csv", index_col=0) #gene list
bar_list = barcodes.index.tolist()
gene_list = genes.index.tolist()
b = adata.layers['spliced']
pancreas_data = pd.DataFrame(data=b.toarray(), index=bar_list, columns=gene_list)
pancreas_data.to_csv(large_root + "/dentategyrus.csv")

#filter and normalize and select genes in this step -- here add ddgs
#scv.pp.filter_and_normalize(adata, **params)

scv.pp.filter_genes(adata, min_shared_counts=20)
scv.pp.normalize_per_cell(adata)
scv.pp.filter_genes_dispersion(adata, n_top_genes=2000)
scv.pp.log1p(adata)

hvglist = adata.var.highly_variable
hvglist.to_csv(large_root + "/dentategyrus_hvg_list.csv")

#calc mean and variance
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

#solve velocity model
#scv.tl.velocity(adata, mode='stochastic', **params)

scv.tl.velocity(adata, mode='dynamical')
scv.tl.velocity_graph(adata)

"""The velocities are projected into a lower-dimensional embedding by translating them into likely cell transitions.
 That is, for each velocity vector we find the likely cell transitions that are in accordance with that direction.
  The probabilities of one cell transitioning into another cell are computed using cosine correlation 
  (between the potential cell transition and the velocity vector) and are stored in a matrix denoted as velocity graph:"""

#adata.write(large_root + 'hvg_model_pancreas.h5ad', compression='gzip')
#adata = scv.read('data/pancreas.h5ad')
scv.pl.velocity_embedding_stream(adata, basis='umap', save='hvg_umap_vel.png')

scv.tl.recover_dynamics(adata)
scv.tl.latent_time(adata)
scv.pl.scatter(adata, color='latent_time', color_map='gnuplot', size=80, save='hvg_latenttime.png')

top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index[:300]
print(top_genes)
topddg = top_genes.tolist()
topddg_frame = pd.DataFrame(topddg)
topddg_frame.to_csv(large_root + "/dentategyrus_hvg_top300latenttime.csv")
scv.pl.heatmap(adata, var_names=top_genes, sortby='latent_time', col_color='clusters', n_convolve=100, save='hvg_topgenes.png')