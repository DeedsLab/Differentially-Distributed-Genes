import pandas as pd
import numpy as np
import os as os
import matplotlib.pyplot as plt
from venn import venn

from matplotlib_venn import venn2, venn3
import matplotlib.pyplot as plt
#large_root = r"/Users/breannesparta/Desktop/results/ddgs/zheng9/marker_genes"

large_root = r"/Users/breannesparta/Desktop/results/ddgs/pancreas"

allmarkers = pd.read_csv(r"/Users/breannesparta/Desktop/results/ddgs/zheng9/marker_genes/PanglaoDB_allcells.tsv", sep='\t', usecols=[0,1,2,3])
cell_type_list = ['Acinar cells','Alpha cells','Beta cells', 'Delta cells', 'Ductal cells', 'Epsilon cells','Gamma (PP) cells', 'Pancreatic progenitor cells', 'Pancreatic stellate cells', 'Peri-islet Schwann cells']
zhengmarkers = allmarkers.loc[allmarkers['cell type'].isin(cell_type_list)]
zheng_genenames = zhengmarkers['official gene symbol'].to_list()
marker_set = set(zheng_genenames)
marker_frame = pd.DataFrame(marker_set)
marker_frame.to_csv("/Users/breannesparta/Desktop/results/ddgs/pancreas/pancreas_Panglao_markergenes.csv")

total_data = pd.read_csv(large_root + "/pancreas.csv", index_col=0, nrows=2) #cell x gene
hvgs = pd.read_csv(large_root + "/pancreas_hvg_list.csv", index_col=0)
ddgs = pd.read_csv(large_root + "/pancreas_ddg_list.csv", index_col=1)
latenttime_hvgs = pd.read_csv(large_root + "/pancreas_hvg_top300latenttime.csv", index_col=1)
latenttime_ddgs = pd.read_csv(large_root + "/pancreas_ddg_top300latenttime.csv", index_col=1)
markers = pd.read_csv(large_root + "/pancreas_Panglao_markergenes.csv", index_col=1,) #gene x cells

all_list = total_data.columns.tolist()
hvg_list = hvgs.index.tolist()
ddg_list = ddgs.index.tolist()
latenttime_hvgs_list = latenttime_hvgs.index.tolist()
latenttime_ddgs_list = latenttime_ddgs.index.tolist()
marker_first = markers.index.tolist()
formated_markers = [str.capitalize(x) for x in marker_first]

marker_overlap = [x for x in formated_markers if x in all_list]

#one marker gene detected in pancreas data; none are in feature sets
#calc jaccard index
def WC_similarity(list1):
    s1 = set(list1)
    s2 = set(gene_list[0]) #this is marker list in this case
    return float(len(s1.intersection(s2)) / len(s1.union(s2)))


title_list = ["marker genes","all_genes","HVGs","top_time_hvgs","DDGs","top_time_ddgs"]
gene_list = [formated_markers, all_list, hvg_list, latenttime_hvgs_list,ddg_list, latenttime_ddgs_list]

# Ven diagram set overlaps
fig = plt.figure(figsize=(15, 8))
ax = fig.add_subplot(111)
venn2([set(gene_list[2]), set(gene_list[4])], set_labels=('HVGs', 'DDGs'))
ax.set_rasterized(True)
plt.savefig(large_root + "/gene_overlap_HVG_DDG.png")
plt.savefig(large_root + "/gene_overlap_HVG_DDG.eps")

# Ven diagram set overlaps
fig = plt.figure(figsize=(15, 8))
ax = fig.add_subplot(111)
venn2([set(gene_list[3]), set(gene_list[5])], set_labels=('top_HVGs', 'top_DDGs'))
ax.set_rasterized(True)
plt.savefig(large_root + "/latenttopgene_overlap_HVG_DDG.png")
plt.savefig(large_root + "/latenttopgene_overlap_HVG_DDG.eps")

#calc jaccard index with WC group
jindex = [WC_similarity(list1) for list1 in gene_list] #marker set in this # bug in j index?? whyyyy
fig = plt.figure(figsize=(12, 8))
ax1 = fig.add_subplot(111)
jindex_frame = pd.DataFrame({'gene_set': title_list, 'AJI': jindex})
jindex_frame.to_csv(large_root + "/pancreas_jaccardindex_w_markergenes.csv")
fs = 12
ax1 = jindex_frame.plot.bar(x='gene_set', y='AJI', rot=0, width = 0.75)
ax1.set_ylabel("Jaccard Index")
ax1.set_xlabel("set compared to marker genes")
ax1.set_title("proportional overlap with marker genes")
ax1.set_ylim(0,1)
plt.grid()
plt.savefig(large_root + "/pancreas_jaccardI_w_markers.eps")