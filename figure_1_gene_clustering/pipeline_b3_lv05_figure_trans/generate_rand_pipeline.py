import pandas as pd
import numpy as np
import scanpy as sc
from sklearn.metrics.cluster import adjusted_rand_score
import matplotlib.pyplot as plt
from sklearn import manifold
from sklearn.neighbors import NearestNeighbors
from sklearn.decomposition import PCA
import smtplib
import os


# read the data from 10x .mtx:
def neighbors(data, k=20):
    # for a given dataset, finds the k nearest neighbors for each point
    nbrs = NearestNeighbors(n_neighbors=k, algorithm='ball_tree').fit(data)
    distances, indices = nbrs.kneighbors(data)
    return indices

def jaccard(A,B):
    # for two sets A and B, finds the Jaccard distance J between A and B
    A = set(A)
    B = set(B)
    union = list(A|B)
    intersection = list(A & B)
    J = ((len(union) - len(intersection))/(len(union)))
    return(J)

def generate_4_rand_samples(root,raw_data_path):
    print("Reading data...")
    large_root=root+"/subclusters/"
    #and "_as_basis" not in item  and "random" not in item])
    new_root=root+"/random_subsets"
    try:
        print("making directory for random subcluster")
        os.makedirs(new_root)
    except OSError:
        print("Creation of the directory %s failed" % new_root)
    total_data_path = raw_data_path
    print("reading total csv files")
    total_data = pd.read_csv(total_data_path)
    total_data =total_data.iloc[:,1:]
    cluster_size=20
    amtclusters = (len([name for name in os.listdir(large_root) if os.path.isfile(os.path.join(large_root,name)) and "asbasisonly" not in name]))
    print(" has "+str(amtclusters) + " subclusters")
    print("Finding High dimensional neighborhoods for all files")
    frame_list=[]
    neighborhoods=[]
    sub_shape = []
    col_labels=['Embedded Dimension']
    for j in range(0,amtclusters):
        if j==0:
            inner_frame=[0]*amtclusters
            inner_n =[0]*amtclusters
        embedded_data= pd.read_csv(large_root+"subcluster_"+str(j)+".csv")
        embedded_data= embedded_data.iloc[:,1:]
        num_cells = embedded_data.shape[0]
        inner_new_root = new_root
        try:
            os.makedirs(inner_new_root)
        except OSError:
            print ("Creation of the directory %s failed" % new_root)
        print("creating " + " random subsample from total data of " +str(num_cells) + " size")
        sub_sample = total_data.sample(n=num_cells)
        sub_sample_2 = total_data.sample(n=num_cells)
        sub_sample_3 = total_data.sample(n=num_cells)
        sub_sample_4 = total_data.sample(n=num_cells)
        sub_sample.to_csv(inner_new_root+"/"+"_sized_1random_subsample_"+str(j)+".csv")
        sub_sample_2.to_csv(inner_new_root+"/"+"_sized_2random_subsample_"+str(j)+".csv")
        sub_sample_3.to_csv(inner_new_root+"/"+"_sized_3random_subsample_"+str(j)+".csv")
        sub_sample_4.to_csv(inner_new_root+"/"+"_sized_4random_subsample_"+str(j)+".csv")
        print("Finding High D Neighborhood for subsample" +str(j))


    ### Write the result to .csv



### if you want to read a loom file:
# adata = sc.read_loom(filename)

# Do a nearest neighbor search:
