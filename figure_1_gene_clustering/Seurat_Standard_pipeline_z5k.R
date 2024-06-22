rm(list=ls())
library(data.table)
library(mclust)
library(fossil)
library(Seurat)
library(dplyr)
library(Matrix)
library(data.table)

siebert_data<-read.table("/home/breanne/sourcedata/Zheng9/zheng5k_dropzeros.csv", sep=",", header = TRUE, row.names=1)
print(dim(siebert_data))
siebert_obj<-CreateSeuratObject(counts=siebert_data, project = "zheng5k")
siebert_obj<-NormalizeData(siebert_obj, normalization.method = "LogNormalize", scale.factor = 10000)
siebert_obj <- FindVariableFeatures(siebert_obj)#Gives 2000 Genes based upon Dispersion See documentation for details
siebert_obj <-ScaleData(siebert_obj)
siebert_obj<-RunPCA(siebert_obj,features = VariableFeatures(siebert_obj), npcs = 50)# This PCA to determine number of PCS in ElbowPlot
ElbowPlot(siebert_obj, ndims=50)
png("/home/breanne/sourcedata/Zheng9/z5k_PCA.png")
dev.off()
siebert_obj<-RunPCA(siebert_obj,features = VariableFeatures(siebert_obj), npcs = 35 )# This is the final value for PCA determination
siebert_obj<-FindNeighbors(siebert_obj,reduction = "pca",graph.name="pca.level", force.recalc = TRUE, dims = 1:35)
siebert_obj<-FindClusters(siebert_obj,graph.name = "pca.level")
reduct<- siebert_obj[["pca"]]
output3<-as.data.frame(as.matrix(Embeddings(reduct)))
PCA_level<-siebert_obj@active.ident#as.numeric(Idents(siebert_obj))# removes vector of cell_barcode and assigned cluster_number
#ordered_PCA_level<-PCA_level[order(row.names(PCA_level)),]
cluster_frame<-data.frame("Cluster_Label"=PCA_level)
fwrite(output3, file ="/home/breanne/sourcedata/Zheng9/zheng5k_PCA_hvgs.csv", row.names = TRUE)#
fwrite(cluster_frame, file ="/home/breanne/sourcedata/Zheng9/zheng5k_Rclusters.csv", row.names = TRUE)#prints vector to csv
siebert_obj <- RunTSNE(siebert_obj, reduction = "pca")
reduct2<- siebert_obj[["tsne"]]
output4<-as.data.frame(as.matrix(Embeddings(reduct2)))
fwrite(output4, file ="/home/breanne/sourcedata/Zheng9/zheng5k_TSNEcoords.csv", row.names = TRUE)
