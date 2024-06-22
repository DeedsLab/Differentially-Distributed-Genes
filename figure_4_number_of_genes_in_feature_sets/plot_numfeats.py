import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.colors
import scipy.stats as ss
import seaborn as sns
from kneed import KneeLocator
from scipy.stats import variation
import scipy.stats as st


large_root = r"/media/timothyhamilton/My Passport1/Tim_Hamilton"
print("Reading data...")
label='feat_sets'
#all data
cytT = pd.read_csv(large_root + r"/Zheng_FACS_pbmcs/cytotoxic_t.csv", index_col=0, usecols = [0,1])
zheng = pd.read_csv(large_root + r"/Zheng_FACS_pbmcs/zheng9_5k/zheng5k_dropzeros.csv", index_col=0,usecols = [0,1])
hydra = pd.read_csv(large_root + r"/R01_Ref_Datasets/Hydra_Results/Hydra_Raw_Cell_Types.csv", nrows=1)
hydra = hydra.T
bladder = pd.read_csv(large_root + r"/Mouse_bladder/mouse_bladder.csv", index_col=0,usecols = [0,1]) #gene by cell
kidney = pd.read_csv(large_root + r"/Mouse_kidney/mkidney.csv", index_col=0,usecols = [0,1])
planaria = pd.read_csv(large_root + r"/Planarian_atlas/planarian.csv", index_col=0,usecols = [0,1]) #nrows=1
#ddgs
#DDG_cytT = pd.read_csv(large_root + r"/Zheng_FACS_pbmcs/cytotoxic_M3Drop.csv", index_col=0, usecols = [0,1])
#DDG_zheng = pd.read_csv(large_root + r"/Zheng_FACS_pbmcs/zheng9_5k/zheng9_5k_DDGs.csv", nrows=1)
#DDG_hydra = pd.read_csv(large_root + r"/Hydra/Hydra_DDG_data.csv", index_col=0,usecols = [0,1])
#DDG_bladder = pd.read_csv(large_root + r"/Mouse_bladder/bladder_DDG_data.csv", index_col=0,usecols = [0,1]) #gene by cell
#DDG_kidney = pd.read_csv(large_root + r"/Mouse_kidney/mkidney.csv", index_col=0,usecols = [0,1])
#m3drop
M3_cytT = pd.read_csv(large_root + r"/Zheng_FACS_pbmcs/cytotoxic_T_M3Drop.csv", nrows=1)
M3_zheng = pd.read_csv(large_root + r"/Zheng_FACS_pbmcs/zheng9_5k/zheng5k_dropzeros_M3Drop.csv", nrows=1)
M3_hydra = pd.read_csv(large_root + r"/R01_Ref_Datasets/Hydra_Results/Hydra_M3Drop_Genes.csv", nrows=1)
M3_bladder = pd.read_csv(large_root + r"/Mouse_bladder/mouse_bladder_M3Drop.csv", nrows=1) #gene by cell
M3_kidney = pd.read_csv(large_root + r"/Mouse_kidney/mkidney_M3Drop.csv", nrows=1)
M3_planaria = pd.read_csv(large_root + r"/Planarian_atlas/planarian_M3Drop.csv", nrows=1) #nrows=1

#NBdrop
NB_cytT = pd.read_csv(large_root + r"/Zheng_FACS_pbmcs/cytotoxic_T_NBUMI.csv", nrows=1)
NB_zheng = pd.read_csv(large_root + r"/Zheng_FACS_pbmcs/zheng9_5k/zheng5k_dropzeros_NBUMI.csv", nrows=1)
NB_hydra = pd.read_csv(large_root + r"/R01_Ref_Datasets/Hydra_Results/Hydra_Raw_Cell_Types_Townes.csv", nrows=1)
NB_bladder = pd.read_csv(large_root + r"/Mouse_bladder/mouse_bladder_NBUMI.csv", nrows=1) #gene by cell
NB_kidney = pd.read_csv(large_root + r"/Mouse_kidney/mkidney_NBUMI.csv", nrows=1)
NB_planaria = pd.read_csv(large_root + r"/Planarian_atlas/planarian_NBUMI.csv", nrows=1) #nrows=1

label = 'feature_sets'
all_list = [cytT,zheng,bladder,kidney,hydra,planaria]
M3_list = [M3_cytT,M3_zheng,M3_bladder,M3_kidney,M3_hydra,M3_planaria]
NB_list = [NB_cytT,NB_zheng,NB_bladder,NB_kidney,NB_hydra,NB_planaria]
label_list = ['cytT','zheng','bladder','kidney','hydra','planaria']

M3_frac = []
NB_frac = []
DDG_list = [3.8,15.2,27.4, 39.1,54.3,55.4]
for data, cname, M3, NB, in zip(all_list,label_list, M3_list, NB_list):
    totalsize = data.shape[0]
    NBsize = (NB.shape[1])/(totalsize)
    M3size = (M3.shape[1])/(totalsize)
    NB_frac.append(NBsize)
    M3_frac.append(M3size)
samp_data = pd.DataFrame({'label_list': label_list, 'M3_data': M3_frac, 'NB_data': NB_frac, 'DDG': DDG_list}, index=label_list)
samp_data.to_csv(large_root + "/Breanne/" + label + "agg_feat_data.csv")

fs = 12
fig = plt.figure(figsize=(5, 8))
ax = fig.add_subplot(111)
ax = samp_data.plot.bar(x='label_list', y='M3_data', rot=0, width=0.75)
ax.set_title("Fraction significant M3Drop genes per sample", fontsize=fs)
ax.set_ylabel("Fraction of genes")
ax.set_xlabel("sample")
ax.set_ylim(0, 0.6)
plt.grid()
plt.savefig(large_root + "/Breanne/" + label + "num_M3drop.png")
plt.savefig(large_root + "/Breanne/" + label + "num_M3Drop.eps")

fs = 12
fig = plt.figure(figsize=(5, 8))
ax = fig.add_subplot(111)
ax = samp_data.plot.bar(x='label_list', y='NB_data', rot=0, width=0.75)
ax.set_title("Fraction significant NBdrop genes per sample", fontsize=fs)
ax.set_ylabel("Fraction of genes")
ax.set_xlabel("sample")
ax.set_ylim(0, 0.6)
plt.grid()
plt.savefig(large_root + "/Breanne/" + label + "num_NBdrop.png")
plt.savefig(large_root + "/Breanne/" + label + "num_NBdrop.eps")

fs = 12
fig = plt.figure(figsize=(5, 8))
ax = fig.add_subplot(111)
ax = samp_data.plot.bar(x='label_list', y='DDG', rot=0, width=0.75)
ax.set_title("Fraction significant DDGs genes per sample", fontsize=fs)
ax.set_ylabel("Fraction of genes")
ax.set_xlabel("sample")
ax.set_ylim(0, 60)
plt.grid()
plt.savefig(large_root + "/" + label + "num_DDGs.png")
plt.savefig(large_root + "/" + label + "num_DDGs.eps")
