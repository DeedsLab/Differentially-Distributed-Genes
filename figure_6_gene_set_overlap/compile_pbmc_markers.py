import pandas as pd
import numpy as np
import os as os
import matplotlib.pyplot as plt

large_root = r"/Users/breannesparta/Desktop/results/ddgs/zheng9/marker_genes"

allmarkers = pd.read_csv(large_root + r"/PanglaoDB_allcells.tsv", sep='\t', usecols=[0,1,2,3])
cell_type_list = ['B cells','Monocytes','NK cells', 'T cells']
zhengmarkers = allmarkers.loc[allmarkers['cell type'].isin(cell_type_list)]
zheng_genenames = zhengmarkers['official gene symbol'].to_list()
marker_set = set(zheng_genenames)
marker_frame = pd.DataFrame(marker_set)
marker_frame.to_csv(large_root + "/zheng_Panglao_markergenes.csv")

markers = pd.read_csv(large_root + "/zheng_Panglao_markergenes.csv", index_col=1,) #gene x cells