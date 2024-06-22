library(M3Drop)
library(dplyr)
library(data.table)
sub_folder<-"/media/timothyhamilton/data1/Tim_Hamilton/Zheng_FACS_pbmcs/zheng9_5k/zheng5k_dropzeros.csv"

datafile<- read.csv(sub_folder, row.names = 1)
print(head(datafile[1,1:6]))

col_list<- colnames(datafile)
#col_1<- as.vector(col_list[2])
#print(col_1)
#tidy_sub<- datafile[,-c(1,2)]
data_mat<- as.data.frame(datafile)
#print(sum(is.na(data_mat)))
#print(colnames(datafile)[1])
norm <- NBumiConvertData(data_mat, is.counts=TRUE)
norm[is.na(norm)]=0
#norm<-t(norm)
model<-NBumiFitModel(norm)
M3Drop_genes <- NBumiFeatureSelectionHighVar(model)
fwrite(as.matrix(M3Drop_genes), file ="/media/timothyhamilton/data1/Tim_Hamilton/Zheng_FACS_pbmcs/zheng9_5k/zheng5k_dropzeros_NBUMI_High_Var_Results.csv", row.names = TRUE)#
print(dim(M3Drop_genes))
gene_list<- names(M3Drop_genes[1:2000])
#subset_list<-append(col_1,gene_list)
gene_data<- as.data.frame(t(as.data.frame(datafile[gene_list,])))
rownames(gene_data)<- col_list
fwrite(gene_data, file ="/media/timothyhamilton/data1/Tim_Hamilton/Zheng_FACS_pbmcs/zheng9_5k/zheng5k_dropzeros_NBUMI_High_Var.csv", row.names = TRUE)#
print(dim(gene_data))
print(head(gene_data)[1,1:6])
