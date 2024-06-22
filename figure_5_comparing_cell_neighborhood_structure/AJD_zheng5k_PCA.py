import pandas as pd
import numpy as np
import numpy.random as random
import matplotlib.pyplot as plt
from sklearn.neighbors import NearestNeighbors
from sklearn.decomposition import PCA

# call functions
def neighbors(data, k=20):
    # for a given dataset, finds the k nearest neighbors for each point
    nbrs = NearestNeighbors(n_neighbors=k + 1, algorithm='ball_tree').fit(data)
    distances, indices = nbrs.kneighbors(data)
    return indices[:, 1:]

def jaccard(A, B):
    # for two sets A and B, finds the Jaccard distance J between A and B
    A = set(A)
    B = set(B)
    union = list(A | B)
    intersection = list(A & B)
    J = ((len(union) - len(intersection)) / (len(union)))
    return (J)

large_root = "/home/breanne/figures/"
label = 'zheng5k'
raw_data = pd.read_csv('/home/breanne/sourcedata/Zheng9/zheng9_5k.csv', index_col = 0)
#sig_data = pd.read_csv(r"C:\Users\bsparta\Desktop\data\R01\sig_data\Duo_8mix_sig_data.csv")
#HVG_list = pd.read_csv(r"C:\Users\bsparta\Desktop\data\R01\Seurat_results\Duo_8mix_Var_features.csv")

#calc AJD as function of PCA embedding for all cells DDG,HVG,all genes
#################################################################33
#clustname = "zheng5k_DDGs"
#total_data = sig_data

#clustname = "zheng5k_HVGs"
#total_data = HVG_list

clustname = "zheng5k"
total_data = raw_data
#################
#calc for each cell type as f(PCA)

cellgroup1 = total_data.filter(regex='^b',axis=1)
cellgroup2 = total_data.filter(regex='^naive_cy',axis=1) #cytotoxic
cellgroup3 = total_data.filter(regex='^mono',axis=1) #monocytes
cellgroup4 = total_data.filter(regex='^regulatory',axis=1) #t
cellgroup5 = total_data.filter(regex='^cd4',axis=1) #thelper
cellgroup6 = total_data.filter(regex='^cd56',axis=1) #nk
cellgroup7 = total_data.filter(regex='^mem',axis=1) #memory t
cellgroup8 = total_data.filter(regex='^naive_t',axis=1) #naive t
cellgroup9 = total_data.filter(regex='^cytotoxic',axis=1) #naive t

data_list = [cellgroup1, cellgroup2, cellgroup3, cellgroup4, cellgroup5, cellgroup6, cellgroup7, cellgroup8, cellgroup9, total_data]
cellname = ['b_cell', 'naive_cytotox', 'monocyte', 'regulatory_T', 'helper_T','NK','memory_T','naive_T','cytotoxic_T','zheng5k']

for cell_group, clabel in zip(data_list,cellname):
    total_data = cell_group.values
    shuffled_gene_copy = total_data.copy()
    for i in range(shuffled_gene_copy.shape[1]):
        random.shuffle(shuffled_gene_copy[:, i])

    mixed_genes = shuffled_gene_copy.copy()
    for i in range(mixed_genes.shape[0]):
        random.shuffle(mixed_genes[i, :])

    min_size = total_data.shape[0]
    if min_size >= 5000:
        step = 100
    if min_size >= 3500:
        step = 50
    elif min_size >= 500:
        step = 20
    elif min_size >= 20:
        step = 5
    else:
        step = 1

    clust_list = [total_data, shuffled_gene_copy, mixed_genes]
    label = ["gene features", "shuffled across cells", "shuffled across genes"]
    dimlist = list(range(1, min_size, step))
    data_list = [dimlist]
    cluster_size = 20
    geneload_array = list()
    gshuffload_array = list()
    cshuffload_array = list()
    for iteration, t in enumerate(clust_list):
        neighborhoods = neighbors(t, k=cluster_size)
        list_array = []
        for dim in dimlist:
            try:
                embedding = PCA(n_components=dim, svd_solver='arpack').fit_transform(t)
                low_D_neighborhood = neighbors(embedding, k=cluster_size)
                jaccard_distances = 0
                for i in range(0, embedding.shape[0], 1):
                    jaccard_distances += jaccard(low_D_neighborhood[i, :], neighborhoods[i, :])
                trial = jaccard_distances / (embedding.shape[0])
            except:
                trial = -1
            list_array.append(trial)
            if iteration == 0:
                geneload_array.append(embedding)
            if iteration == 1:
                gshuffload_array.append(embedding)
            if iteration == 2:
                cshuffload_array.append(embedding)
        data_list.append(list_array)
    nparray = np.asarray(data_list)
    nparray = np.transpose(nparray)
    col_labels = ['Embedded Dimension'] + label
    frame_of_data = pd.DataFrame(nparray, columns=col_labels)
    frame_of_data.to_csv(large_root + "/" + clustname + clabel + "AJD.csv")

    fig = plt.figure(figsize=(15, 8))
    ax1 = fig.add_subplot(111)
    ax1.set_title(
        "eigengenes: " + clustname + clabel)
    ax1.set_ylabel("Average Jaccard Distance")
    ax1.set_ylim((0, 1))
    ax1.set_xlim((0, embedding.shape[1]))
    headers = frame_of_data.keys()
    label = headers[2]
    frame_of_data.plot(x=headers[0], y=headers[1:], ax=ax1, legend=True)
    plt.savefig(large_root + "/" + clustname + clabel + "AJDgenes.png")
    plt.savefig(large_root + "/" + clustname + clabel + "AJDgenes.eps")

    np.savez(large_root + "/" + clustname + clabel + "gene clusters gene features", *geneload_array)
    #np.savez(large_root + "/" + clustname + clabel + "gene clusters shuffled across gene", *gshuffload_array)
    #np.savez(large_root + "/" + clustname + clabel + "gene clusters shuffled across cells", *cshuffload_array)