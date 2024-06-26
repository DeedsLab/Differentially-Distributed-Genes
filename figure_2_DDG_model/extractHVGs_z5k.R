rm(list=ls())
library(data.table)
library(mclust)
library(fossil)
library(Seurat)
library(dplyr)
library(Matrix)
library(data.table)

#siebert_data<-read.table("/home/breanne/sourcedata/Zheng9/zheng9.csv", sep=",", header = TRUE, row.names=1)
siebert_data<-read.table("/home/breanne/sourcedata/Zheng9/zheng5k_dropzeros.csv", sep=",", header = TRUE, row.names=1)
siebert_data_2<- t(siebert_data)
print(dim(siebert_data))
siebert_obj<-CreateSeuratObject(counts=siebert_data_2, project = "zheng5k")
siebert_obj<-NormalizeData(siebert_obj, normalization.method = "LogNormalize", scale.factor = 10000)
#siebert_obj <-ScaleData(siebert_obj)
siebert_obj <- FindVariableFeatures(siebert_obj)#Gives 2000 Genes based upon Dispersion See documentation for details
feat<- VariableFeatures(siebert_obj)
siebert_data_var<- siebert_data[feat,]
fwrite(siebert_data_var, file ="/home/breanne/sourcedata/Zheng9/zheng5k_HVGs.csv", row.names = TRUE)
