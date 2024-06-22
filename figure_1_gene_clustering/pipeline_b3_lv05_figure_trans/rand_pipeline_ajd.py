import pandas as pd
import numpy as np
import scanpy as sc
from sklearn.metrics.cluster import adjusted_rand_score
import matplotlib.pyplot as plt
import os
import sys
from sklearn import manifold
from sklearn.neighbors import NearestNeighbors
from sklearn.decomposition import PCA
import smtplib
import os, os.path
pathlist=["corl279","dms454", "dms53", "h1048","h524","h69", "h82", "h841", "mixbackup", "mix"]
#merge the subclusters
def neighbors (data,k=20):
	nbrs=NearestNeighbors(n_neighbors=k+1, algorithm = 'ball_tree').fit(data)
	distances,indices = nbrs.kneighbors(data)
	return indices

def jaccard (A,B):
	A=set(A)
	B=set(B)
	union=list(A|B) #separated sets by a pipeline
	intersection = list(A & B)
	J = ((len(union)-len(intersection))/(len(union)))
	return J

def rand_ajd_generate(cluster_size,root):
	rootpath = root+"/random_subsets/"
	print("well, good luck, :)")
	clusters= ([name for name in os.listdir(rootpath) if os.path.isfile(os.path.join(rootpath,name))])
	print("there are "+ str(len(clusters)) + " clusters in set")
	frame_list=[]
	neighborhoods=[]
	sub_shape=[]
	nameit=[]
	for j in clusters:
		print("performing operations to create"+"subcluster_"+str(j)+"asbasisonly.csv")
		embedded_data = pd.read_csv(rootpath+str(j))
		embedded_data=embedded_data.iloc[:,1:]
		print("Finding High_D Neighborhood")
		high_D_neighborhood= neighbors(embedded_data,k=cluster_size)
		frame_list.append(embedded_data)
		neighborhoods.append(high_D_neighborhood)
		sub_shape.append(embedded_data.shape[0])
		nameit.append(rootpath+str(j)+"basis.csv")
	for b,c,d in zip(neighborhoods,frame_list,nameit):
		basis_data=c
		data_list=[]
		min_size= basis_data.shape[0]
		if min_size>=3500:
			step = 50
		elif min_size>=500:
			step = 20
		elif min_size>=20:
			step = 5
		else:
			step = 1
		for dim in range(2,min_size,step):
			list_array=[dim]
			try:
				print("Generating Embedding...")
				embedding= PCA(n_components=dim,svd_solver='auto').fit_transform(basis_data)
				print("Finding Low_D Neighborhood...")
				low_D_neighborhood= neighbors(embedding,k=cluster_size)
				print("Calculating Jaccard Distances...")
				jaccard_distances=0
				for s in range(0,embedding.shape[0],1):
					jaccard_distances+=jaccard(low_D_neighborhood[s,:],b[s,:]) #ask tim about this
				trial = jaccard_distances/(embedding.shape[0])
			except:
				print("Did not work with latent dim:" + str(dim) + "-1 put in as placeholder and will remove when plotting")
				trial= -1
			list_array.append(trial)
			data_list.append(list_array)
			nparray=np.asarray(data_list)
			col_labels=['Embedded Dimension','PCA']
			frame_of_data = pd.DataFrame(nparray, columns=col_labels)
			print("Making .csv")
				#dot=a.index(".") #ask tim about this
			frame_of_data.to_csv(d)
			print("Done with one Dimension")

