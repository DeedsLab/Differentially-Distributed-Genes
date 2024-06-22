import pandas as pd
import numpy as np
import os as os
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy.stats import variation
from scipy.stats import binom
import numpy.random as random
import math

#select 2000 random genes

large_root = '/home/breannesparta/gcc_data/mTandB'
rep_list = [1,2,3,4,5,6,7,8,9,10]
#rep_list = [1]

#cell x gene
data = pd.read_csv(large_root + r"/BandT.csv", index_col=0)
cell_headers = data["Cell_Type"]
data = data.drop(columns=["Cell_Type"])

for rep in rep_list:
    rep = str(rep)
    ngenes = data.shape[1]
    print(data.shape)

    randnums= np.random.randint(1,ngenes,2000)
    bcells = data.iloc[:,randnums]

    bcells.insert(0, "Cell_Type", cell_headers)
    bcells.to_csv("/home/breannesparta/ddgs/BandT/random2k/2k_rand_BandT_rep_" + rep + ".csv")