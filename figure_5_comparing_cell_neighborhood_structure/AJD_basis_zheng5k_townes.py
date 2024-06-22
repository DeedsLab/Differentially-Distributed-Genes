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

large_root = r"/media/timothyhamilton/My Passport1/Tim_Hamilton/Zheng_FACS_pbmcs/zheng9_5k"
label = 'zheng9_5k'

total_data = pd.read_csv(r"/media/timothyhamilton/My Passport1/Tim_Hamilton/Zheng_FACS_pbmcs/zheng9_5k/zheng5k_dropzeros.csv", index_col=0) #genes x cells
WC = pd.read_csv(r"/media/timothyhamilton/My Passport1/Tim_Hamilton/Zheng_FACS_pbmcs/zheng9_5k/zheng9_5k_WC.csv", index_col=0) #g x c
HVG = pd.read_csv(r"/media/timothyhamilton/My Passport1/Tim_Hamilton/Zheng_FACS_pbmcs/zheng9_5k/zheng5k_HVGs.csv", index_col=0) #g x c
HVG = HVG.replace(np.nan, 0)
DDG = pd.read_csv(r"/media/timothyhamilton/My Passport1/Tim_Hamilton/Zheng_FACS_pbmcs/zheng9_5k/zheng9_5k_DDGs.csv", index_col=0) #g x c
townes = pd.read_csv(r"/media/timothyhamilton/My Passport1/Tim_Hamilton/Zheng_FACS_pbmcs/zheng9_5k/zheng5k_dropzeros_SCRY.csv", index_col=0)#c x g
mdrop = pd.read_csv(r"/media/timothyhamilton/My Passport1/Tim_Hamilton/Zheng_FACS_pbmcs/zheng9_5k/zheng5k_dropzeros_M3Drop.csv", index_col=0) # c x g
nbdrop = pd.read_csv(r"/media/timothyhamilton/My Passport1/Tim_Hamilton/Zheng_FACS_pbmcs/zheng9_5k/zheng5k_dropzeros_NBUMI.csv", index_col=0)#c x g
nbdisp = pd.read_csv(r"/media/timothyhamilton/My Passport1/Tim_Hamilton/Zheng_FACS_pbmcs/zheng9_5k/zheng5k_dropzeros_NBUMI_High_Var.csv", index_col=0) # c x g

total_data = total_data.T
WC = WC.T
HVG = HVG.T
DDG = DDG.T

total_data = total_data.values
WC = WC.values
HVG = HVG.values
DDG = DDG.values
townes = townes.values
mdrop = mdrop.values
nbdrop = nbdrop.values
nbdisp = nbdisp.values

cluster_size = 20
total_neighborhood = pd.DataFrame(neighbors(total_data, k=cluster_size))
WC_neighborhood = pd.DataFrame(neighbors(WC, k=cluster_size))
HVG_neighborhood = pd.DataFrame(neighbors(HVG, k=cluster_size))
DDG_neighborhood = pd.DataFrame(neighbors(DDG, k=cluster_size))
townes_neighborhood = pd.DataFrame(neighbors(townes, k=cluster_size))
mdrop_neighborhood = pd.DataFrame(neighbors(mdrop, k=cluster_size))
nbdisp_neighborhood = pd.DataFrame(neighbors(nbdisp, k=cluster_size))
nbdrop_neighborhood = pd.DataFrame(neighbors(nbdrop, k=cluster_size))

neighborhoods = [total_neighborhood, WC_neighborhood, HVG_neighborhood, DDG_neighborhood, townes_neighborhood, mdrop_neighborhood, nbdisp_neighborhood, nbdrop_neighborhood]
data_names = ['all_zheng5k','WC', 'HVGs','DDGs','townes','m3drop','nbdisp','nbdrop']

xlabels = list(data_names)
ylabels = list(data_names)
AJDdata = pd.DataFrame(data_names, columns=["First"]).set_index('First', drop=True)

#for neighbors, glabel in zip(neighborhoods, data_names):
neighbors = townes_neighborhood
glabel = 'townes'
list_array = []
for neighbors2 in neighborhoods:
    jaccard_dist = 0
    for s in range(0, total_neighborhood.shape[0], 1):
        jaccard_dist += jaccard(neighbors.iloc[s, :], neighbors2.iloc[s, :])
        trial = (jaccard_dist/total_neighborhood.shape[0])
    list_array.append(trial)
AJDdata[glabel] = list_array
AJDd = AJDdata
AJDd.to_csv(large_root + "/" + label + "_AJDcompare_basist.csv")
#AJDdata.to_csv(large_root + "/" + label + "_AJDcompare_basis.csv")

fs = 12
fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)
sns.heatmap(AJDdata, vmin = 0, vmax = 1, square=True)
ax.set_title("AJD: bases comparison \n" + label, fontsize=fs)
plt.tight_layout()
plt.savefig(large_root + "/" + label + "_AJDcompare_basist.png")
plt.savefig(large_root + "/" + label + "_AJDcompare_basist.eps")
