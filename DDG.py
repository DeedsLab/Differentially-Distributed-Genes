import pandas as pd
import numpy as np
import time
import math
import scipy.stats
import sys

print("Loading in csv")

large_root = "100_rep_1zheng9_resamp.csv"
capture_probability = 0.05

# read csv or scanpy obj. Note that this script works for genes as columns and cells as rows.
bcells1 = pd.read_csv(large_root, index_col=0, engine = "c")
bcells1 = bcells1.loc[:, (bcells1 != 0).any()]  # drop genes that are not detected in any cell

print("Calculating relevant metadata")

#Create variables for number of genes samples and cells used
cell_count, gene_count = bcells1.shape
genes = bcells1.columns
#Create vectorized representations relevant cellular information on a gene by gene basis
cells_with_mRNAs = bcells1.astype(bool).sum(axis=0)
average_expression = round(bcells1.sum(axis=0)/cell_count, 6)
mRNA_amount = average_expression/capture_probability
no_mRNA_in_cell = np.exp(mRNA_amount * np.log(1 -  capture_probability))
cells_with_no_mRNAs = cell_count - cells_with_mRNAs

#Create an array to store p-values

p_val_array = np.zeros(gene_count)

start_time = time.time()

print("starting P val calculations")

#create an array to store p-values, calculate p values

for n in range(0,gene_count):
    count_array = np.arange(cell_count - cells_with_mRNAs[n], cell_count)
    p_val_array[n] = np.sum(scipy.stats.binom.pmf(count_array, cell_count, no_mRNA_in_cell[n]))



""" The rounding present may impact the significance of the selected genes.


p_val_array = (1000000 * p_val_array).astype(int) / 1000000

p_val_array = np.round(p_val_array, 6)


"""

p_val_array = np.round(p_val_array, 6)

#arrange p values and corresponding genes into new array

final_array = np.vstack((genes,p_val_array)).T

#Save final array
pd.DataFrame(final_array).to_csv(large_root[:-4]+"_PVALS.csv", index = False, header = False)

print("Finished")

sys.exit()








