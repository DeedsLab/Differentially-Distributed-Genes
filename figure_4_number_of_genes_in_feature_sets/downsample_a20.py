import pandas as pd
import numpy as np
import os as os
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy.stats import variation
from scipy.stats import binom
import numpy.random as random
import math


def do_bernoulli_trial(n):
    """Perform n Bernoulli trials with success probability p
    and return number of successes."""
    # Initialize number of successes: n_success
    # print(n)
    n_success = 0
    for i in range(n):
        random_number = np.random.random()
        # If less than p, it's a success  so add one to n_success
        if random_number < problist[0]:
            n_success += 1
    return n_success


large_root = "/home/breannesparta/ddgs/A20_3T3/cellgene"
rep_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

# capture probability
problist = [0.9]

# cell x gene
label = "a20_3t3"
data = pd.read_csv(large_root + "/A20_3T3.csv", index_col=0)
# data = data.T
print(data.head())  # cell by gene data
# gene_list = data.columns.to_list()
data = data.drop(columns=["Cell_Type"])

for rep in rep_list:
    reps = str(rep)
    new_data = data.copy()
    # new_data = new_data.round()
    new_data = new_data.astype('Int64')
    for prob in problist:
        probs = str(prob)
        #downsampled = new_data.applymap(lambda x: do_bernoulli_trial(x))[new_data > 0] #filter > 0 gives bug idk
        downsampled = new_data.applymap(lambda x: do_bernoulli_trial(x))
        downsampled = downsampled.fillna(0)
        print(downsampled.head())
        downsampled.to_csv(large_root + "/" + label + "downsample" + probs + "reps" + reps + ".csv")
