import pandas as pd
import numpy as np
import scanpy as sc
import numpy.random as random
from sklearn.metrics.cluster import adjusted_rand_score
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn import manifold
import scipy.stats as ss

from sklearn.neighbors import NearestNeighbors
from sklearn.decomposition import PCA
from sklearn.decomposition import NMF
from sklearn.decomposition import SparsePCA
from sklearn.decomposition import FastICA
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import pairwise_distances

import smtplib
import os
import sys

# read the data from 10x .mtx:
def neighbors(data, k=20):
    # for a given dataset, finds the k nearest neighbors for each point
    nbrs = NearestNeighbors(n_neighbors=k+1, algorithm='ball_tree').fit(data)
    distances, indices = nbrs.kneighbors(data)
    return indices[:,1:]


def max_radius_used(data):
    dis_mat= pairwise_distances(data)
    return  dis_mat

def make_Linear_Multivar_3_Gauss(n_points, n_genes, sep):
    mean_vec=[100]*n_genes
    cov_mat= np.identity(n_genes)
    array_1 = np.random.multivariate_normal(mean_vec,cov_mat,size=n_points)
    label_list=["Population_1"]*array_1.shape[0]
    for multi in range(1,3):
        new_mean_vec= np.asarray([100+multi*sep]*n_genes)
        new_array= np.random.multivariate_normal(new_mean_vec,cov_mat,size=n_points)
        label_list+=["Population_"+str(multi+1)]*new_array.shape[0]
        array_1= np.append(array_1,new_array,axis=0)
    print(array_1.shape)
    new_frame=pd.DataFrame(array_1, columns= ["Gene_"+str(i) for i in range(array_1.shape[1])])
    new_frame.insert(0,"Cell_Type",label_list)
    return new_frame

def make_Marker_Multivar_3_Gauss(n_points, n_genes, basal_val, marker_val):
    mean_vec=[marker_val]+[basal_val]*(n_genes-1)
    cov_mat= np.identity(n_genes)
    array_1 = np.random.multivariate_normal(mean_vec,cov_mat,size=n_points)
    label_list=["Population_1"]*array_1.shape[0]
    for multi in range(1,3):
        #pos 1 corresponds to pop 1;
        new_mean_vec= [marker_val if i== multi else basal_val for i in range(n_genes)]#np.asarray([100+multi*sep]*n_genes)
        new_array= np.random.multivariate_normal(new_mean_vec,cov_mat,size=n_points)
        label_list+=["Population_"+str(multi+1)]*new_array.shape[0]
        array_1= np.append(array_1,new_array,axis=0)
    print(array_1.shape)
    new_frame=pd.DataFrame(array_1, columns= ["Gene_"+str(i) for i in range(array_1.shape[1])])
    new_frame.insert(0,"Cell_Type",label_list)
    return new_frame

label = 'GMM'

raw_mat= make_Marker_Multivar_3_Gauss(4000,200,10,100)
#make each gene an int in the simulation

#save marker gene model

# capture probability
problist = [0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1]

for rep in rep_list:
    reps = str(rep)
    new_data = raw_mat.copy()
    new_data = new_data.round()
    new_data = new_data.astype('Int64')
    for prob in problist:
        probs = str(prob)
        downsampled = new_data.applymap(lambda x: do_bernoulli_trial(x))
        downsampled = downsampled.fillna(0)
        print(downsampled.head())
        downsampled.to_csv(large_root + "/" + label + "downsample" + probs + "reps" + reps + ".csv")

#make distance matrix
"""dist_mat= max_radius_used(raw_mat.iloc[:,1:].values)
print(np.amin(dist_mat[np.nonzero(dist_mat)]))
print(np.amax(dist_mat[np.nonzero(dist_mat)]))
label_1=[]
label_2=[]
dist_list=[]
for rows in range(dist_mat.shape[0]):
    for cols in range(rows+1, dist_mat.shape[1]):
        label_1.append(rows)
        label_2.append(cols)
        dist_list.append(dist_mat[rows,cols])
        label_1.append(cols)
        label_2.append(rows)
        dist_list.append(dist_mat[rows,cols])
result_Frame= pd.DataFrame({"Col_1": label_1, "Col_2": label_2, "Col_3": dist_list})
result_Frame.to_csv(frame_out_path, sep='\t', header=False, index= False)"""
