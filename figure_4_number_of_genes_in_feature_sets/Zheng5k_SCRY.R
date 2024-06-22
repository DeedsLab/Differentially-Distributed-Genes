library(scry)
library(dplyr)
library(data.table)
sub_folder<-"/media/timothyhamilton/My Passport1/Tim_Hamilton/Zheng_FACS_pbmcs/zheng9_5k/zheng5k_dropzeros.csv"

datafile<- read.csv(sub_folder, row.names = 1)
print(head(datafile[1,1:6]))

col_list<- colnames(datafile)
row_list<- rownames(datafile)
#col_1<- as.vector(col_list[2])
#print(col_1)
#tidy_sub<- datafile[,-c(1,2)]
data_mat<- as.data.frame(datafile)
#print(sum(is.na(data_mat)))
#print(colnames(datafile)[1])
my_genes<-devianceFeatureSelection(as.matrix(data_mat))
M3Drop_genes<-data.frame(row_list,my_genes)
M3Drop_genes<-M3Drop_genes[order(-my_genes),]
fwrite(M3Drop_genes, file ="/media/timothyhamilton/My Passport1/Tim_Hamilton/Zheng_FACS_pbmcs/zheng9_5k/zheng5k_dropzeros_SCRY_Results.csv", row.names = TRUE)#
print(dim(M3Drop_genes))
gene_list<- as.vector(M3Drop_genes[1:2000,1])
#subset_list<-append(col_1,gene_list)
gene_data<- as.data.frame(t(as.data.frame(datafile[gene_list,])))
rownames(gene_data)<- col_list
fwrite(gene_data, file ="/media/timothyhamilton/My Passport1/Tim_Hamilton/Zheng_FACS_pbmcs/zheng9_5k/zheng5k_dropzeros_SCRY.csv", row.names = TRUE)#
print(dim(gene_data))
print(head(gene_data)[1,1:6])
