library(M3Drop)
library(dplyr)
library(data.table)
sub_folder<-"/media/timothyhamilton/My Passport1/Tim_Hamilton/Planarian_atlas/planarian.csv"

datafile<- read.csv(sub_folder, row.names = 1)

datafile<- t(datafile)
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
M3Drop_genes <-  NBumiFeatureSelectionCombinedDrop(model,method="fdr", qval.thres=0.05, suppress.plot=FALSE)
fwrite(M3Drop_genes, file ="/media/timothyhamilton/My Passport1/Tim_Hamilton/Planarian_atlas/planarian_NBUMI_Results.csv", row.names = TRUE)#
print(dim(M3Drop_genes))
gene_list<- as.vector(M3Drop_genes[,1])
#subset_list<-append(col_1,gene_list)
gene_data<- as.data.frame(t(as.data.frame(datafile[gene_list,])))
rownames(gene_data)<- col_list
fwrite(gene_data, file ="/media/timothyhamilton/My Passport1/Tim_Hamilton/Planarian_atlas/planarian_NBUMI.csv", row.names = TRUE)#
print(dim(gene_data))
print(head(gene_data)[1,1:6])
