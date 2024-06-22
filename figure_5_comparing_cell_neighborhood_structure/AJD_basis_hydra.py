import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.colors
import scipy.stats as ss
import seaborn as sns
from kneed import KneeLocator
from scipy.stats import variation
from sklearn.neighbors import NearestNeighbors
from sklearn.decomposition import PCA
import math
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

large_root = "/media/timothyhamilton/data1/Tim_Hamilton/Hydra"
label = 'Hydra'
#AJD calc on cell by gene data
total_data = pd.read_csv("/media/timothyhamilton/data1/Tim_Hamilton/Hydra/Hydra_full_Raw.csv", index_col=0) #cell x genes
HVG = pd.read_csv("/media/timothyhamilton/data1/Tim_Hamilton/Hydra/Hydra_HVGs.csv", index_col=0) #c x g
DDG = pd.read_csv("/media/timothyhamilton/data1/Tim_Hamilton/Hydra/Hydra_DDG_data.csv", index_col=0)# g x c
townes = pd.read_csv("/media/timothyhamilton/data1/Tim_Hamilton/Hydra/Hydra_Raw_Cell_Types_Townes.csv", index_col=0)#c x g
mdrop = pd.read_csv("/media/timothyhamilton/data1/Tim_Hamilton/Hydra/Hydra_M3Drop_Genes.csv", index_col=0)
nbdrop = pd.read_csv("/media/timothyhamilton/data1/Tim_Hamilton/R01_Ref_Datasets/Hydra_Results/Hydra_NBUMIDrop_Genes.csv", index_col=0)
nbdisp = pd.read_csv("/media/timothyhamilton/data1/Tim_Hamilton/R01_Ref_Datasets/Hydra_Results/Hydra_NBUMIDrop_Genes_High_Var.csv", index_col=0)
nbdrop = nbdrop.drop(['Cell_Type'], axis = 1)
nbdisp = nbdisp.drop(['Cell_Type'], axis = 1)
townes = townes.drop(['Cell_Type'], axis = 1)
mdrop = mdrop.drop(['Cell_Type'], axis = 1)
DDG = DDG.T

total_data = total_data.values
HVG = HVG.values
DDG = DDG.values
townes = townes.values
mdrop = mdrop.values
nbdrop = nbdrop.values
nbdisp = nbdisp.values

cluster_size = 20
total_neighborhood = pd.DataFrame(neighbors(total_data, k=cluster_size))
HVG_neighborhood = pd.DataFrame(neighbors(HVG, k=cluster_size))
DDG_neighborhood = pd.DataFrame(neighbors(DDG, k=cluster_size))
townes_neighborhood = pd.DataFrame(neighbors(townes, k=cluster_size))
mdrop_neighborhood = pd.DataFrame(neighbors(mdrop, k=cluster_size))
nbdisp_neighborhood = pd.DataFrame(neighbors(nbdisp, k=cluster_size))
nbdrop_neighborhood = pd.DataFrame(neighbors(nbdrop, k=cluster_size))

neighborhoods = [total_neighborhood, HVG_neighborhood, DDG_neighborhood, townes_neighborhood, mdrop_neighborhood, nbdisp_neighborhood, nbdrop_neighborhood]
data_names = ['all_hydra','HVGs','DDGs','townes','m3drop', 'nbdisp','nbdrop']

xlabels = list(data_names)
ylabels = list(data_names)
AJDdata = pd.DataFrame(data_names, columns=["First"]).set_index('First', drop=True)

print("calculating AJD...")

for neighbors, glabel in zip(neighborhoods, data_names):
    list_array = []
    for neighbors2 in neighborhoods:
        jaccard_dist = 0
        for s in range(0, total_neighborhood.shape[0], 1):
            jaccard_dist += jaccard(neighbors.iloc[s, :], neighbors2.iloc[s, :])
            trial = (jaccard_dist/total_neighborhood.shape[0])
        list_array.append(trial)
    AJDdata[glabel] = list_array
    AJDd = AJDdata
    AJDd.to_csv(large_root + "/" + label + "_AJDcompare_basis.csv")
AJDdata.to_csv(large_root + "/" + label + "_AJDcompare_basis.csv")

fs = 12
fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)
sns.heatmap(AJDdata, vmin = 0, vmax = 1, square=True)
ax.set_title("AJD: bases comparison \n" + label, fontsize=fs)
plt.tight_layout()
plt.savefig(large_root + "/" + label + "_AJDcompare_basis.png")
plt.savefig(large_root + "/" + label + "_AJDcompare_basis.eps")
