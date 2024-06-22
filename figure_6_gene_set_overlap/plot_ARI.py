
#about:
# load lists of genes
# format gene names if necessary
# set reference gene list in Jaccard index function
# calculate proportional overlap across gene lists


#format input as list
#ensure gene names are in same format
#load lists of genes for comparison
#WC =

def rename(A):
    # for gene list A return gene names with spaces instead of - or . characters
    newlist = []
    for gene in A:
        newgene = str(gene)
        a = newgene.replace('.', '-') #R changes - to .
        newlist.append(a)
    all_headers = np.array(newlist)
    return (all_headers)

#rename
WC_set= rename(WC.columns.tolist())

#input names and list of gene sets for jaccard index calculation
title_list = ["WC","DDGs","M3Drop","Townes","NBdrop","NBdisp","HVGs","marker genes"]
gene_list = [WC_set, DDG_set, M3Drop_set, Townes_set,nbdrop_set, nbdisp_set,HVG_set]

## input reference dataset for ARI into jaccard index function
def DEG_similarity(list1):
    s1 = set(list1)
    s2 = set(gene_list[0]) #this is Wilcoxon list in this case
    return float(len(s1.intersection(s2)) / len(s1.union(s2)))

#calc jaccard index with WC group
jindex = [DEG_similarity(list1) for list1 in gene_list]
fig = plt.figure(figsize=(12, 8))
ax1 = fig.add_subplot(111)
jindex_frame = pd.DataFrame({'gene_set': title_list, 'AJI': jindex})
fs = 12
#plot results
ax1 = jindex_frame.plot.bar(x='gene_set', y='AJI', rot=0, width = 0.75)
ax1.set_ylabel("Jaccard Index")
ax1.set_xlabel("set compared to WC genes")
ax1.set_title("proportional overlap with WC genes")
ax1.set_ylim(0,1)
plt.grid()

#plot ven diagram
for gene_group, glabel in zip(gene_list, title_list):
    fig = plt.figure(figsize=(15, 8))
    ax = fig.add_subplot(111)
    venn2([set(gene_list[0]), set(gene_group)], set_labels=('WCs', glabel))
    ax.set_rasterized(True)