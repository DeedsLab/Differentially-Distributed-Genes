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
# non random pathlist=["corl279","dms454", "dms53", "h1048","h524","h69", "h82", "h841", "mixbackup", "mix"]

#merge the subclusters

def plot_ajd_generate(cluster_size,root):
	rootpath1=root
	rootpath2=root+"/random_subsets"
	print("well, good luck, :)")
#		if israndom==False:
#			sourcedata = sorted([name for name in os.listdir(rootpath) if os.path.isfile(os.path.join(rootpath,name)) and "Fixed" in name])
#		elif israndom==True:
#			sourcedata= sorted([name for name in os.listdir(rootpath+i+"/") if os.path.isfile(os.path.join(rootpath+i+"/",name)) and "_sized_random_subsample" not in name])
	size_data=len(sorted([name for name in os.listdir(rootpath1+"/subclusters/") if os.path.isfile(os.path.join(rootpath1+"/subclusters/",name)) and "asbasisonly" in name]))
	for i in range(size_data):
		sourcedata1=sorted([name for name in os.listdir(rootpath1+"/subclusters/") if os.path.isfile(os.path.join(rootpath1+"/subclusters/",name)) and "_"+str(i)+"asbasisonly" in name])
		sourcedata2=sorted([name for name in os.listdir(rootpath2+"/") if os.path.isfile(os.path.join(rootpath2+"/",name)) and "random_subsample_"+str(i)+".csvbasis" in name])
		print(sourcedata1)
		for k in sourcedata1:
			counter=i
			print("creating ajd vs. dimension graph for ..." )
			df1 = pd.read_csv(rootpath1+"/subclusters/"+str(k))
			df2list=[]
			for m in sourcedata2:
				df2 = pd.read_csv(rootpath2+"/"+str(m))
				df2list.append(df2)
			headers= df1.keys()
			print(len(df2list))
			headerslist=[i.keys()[1:] for i in df2list]
			fig=plt.figure(figsize=(15,8))
			ax1=fig.add_subplot(111)
			for heads in headers[2:]:
				print("now")
				ax1.set_title("Cytotrace Random vs. NonRandom AJD Subcluster_" +str(i))
				ax1.set_ylabel("Average Jaccard Distance")
				ax1.set_xlabel("Embedded Dimension")
				df1.plot(x=headers[1],y=heads,ax=ax1,label="Regular")
				count2=0
			for o in headerslist:
				df22=df2list[count2]
				headers2=o
				df22.plot(x="Embedded Dimension",y="PCA",ax=ax1,label="Random"+str(count2))
				count2=count2+1
			plt.savefig(rootpath1+"/"+str(counter)+"_subcluster_both.png")
			plt.close()
			counter=counter+1