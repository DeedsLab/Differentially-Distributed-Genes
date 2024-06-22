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
AJDdata = pd.read_csv("/media/timothyhamilton/data1/Tim_Hamilton/Hydra/HydraAJDcompare_basis.csv", index_col=0)

fs = 12
fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)
sns.heatmap(AJDdata, vmin = 0, vmax = 1, square=True)
ax.set_title("AJD: bases comparison \n" + label, fontsize=fs)
plt.tight_layout()
plt.savefig(large_root + "/" + label + "AJDcompare_basis.png")
plt.savefig(large_root + "/" + label + "AJDcompare_basis.eps")
