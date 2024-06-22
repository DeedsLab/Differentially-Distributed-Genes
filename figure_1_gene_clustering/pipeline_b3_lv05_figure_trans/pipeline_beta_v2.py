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
import generate_rand_pipeline
import function_pipeline_ajd
import rand_pipeline_ajd
import pipeline_ajd_plot
import pipeline_louvain
import multiprocessing as mp
from multiprocessing import Process
from multiprocessing import Pool
import subprocess
import time
print("please enter your root that contains your raw file folder")
c=str(input())
print("please enter your raw_data file name")
raw=str(input())
print("now please enter your raw_data delimiter style")
delim=str(input())
start_time=time.time()
if __name__ == '__main__':
	p0=Process(target=pipeline_louvain.louvain_subclustering,args=(c,raw,delim))
	p1=Process(target=function_pipeline_ajd.ajd_generate,args=(20,c))
	p2=Process(target=generate_rand_pipeline.generate_4_rand_samples,args=(c,raw))
	p3=Process(target=rand_pipeline_ajd.rand_ajd_generate,args=(20,c))
	p4=Process(target=pipeline_ajd_plot.plot_ajd_generate,args=(20,c))
#	p4=Process(target=)
	p0.start()
	p0.join()
	p1.start()
	p2.start()
	p2.join()
	p3.start()
	p3.join()
	p1.join()
	p4.start()
	p4.join()
print(' this is how long it took in seconds' +str(time.time()-start_time))
#	p4.start()
#function_pipeline_ajd.ajd_generate(20,c)
#generate_rand_pipeline.generate_4_rand_samples(c,raw)
#rand_pipeline_ajd.rand_ajd_generate(20,c)
#if __name__ =='__main__':
#	p1=Process(function_pipeline_ajd.ajd_generate(20,c))
#	p2=Process(target=generate_rand_pipeline.generate_4_rand_samples,args=(c,raw))
#	#p3=Process(target=rand_pipeline_ajd.rand_ajd_generate,args=(20,c))
#	p1.start()
#	p2.start()
#	print(p1.is_alive(),p2.is_alive())
#	p1.join()
#	p2.join()
#let's do this