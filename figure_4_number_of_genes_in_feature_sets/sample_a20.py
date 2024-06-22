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

large_root = '/home/breannesparta/gcc_data/A20_3T3'
rep_list = [1,2,3,4,5,6,7,8,9,10]
#rep_list = [1]

#cell x gene
data = pd.read_csv(large_root + r"/A20_3T3.csv", index_col=0)
cell_headers = data["Cell_Type"]
data = data.drop(columns=["Cell_Type"])

for rep in rep_list:
    rep = str(rep)
    ngenes = data.shape[1]

    randnums= np.random.randint(1,ngenes,2000)
    bcells = data.iloc[:,randnums]

    bcells.insert(0, "Cell_Type", cell_headers)
    bcells.to_csv("/home/breannesparta/ddgs/A20_3T3/cellgene/2k_rand_A20_3T3_rep_" + rep + ".csv")