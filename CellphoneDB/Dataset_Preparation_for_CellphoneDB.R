library("biomaRt")

loc_anndata<-"Path_to_csv_output"
X<-as.matrix(read.table(paste0(loc_anndata,"/X.csv"),sep = ","))
var<-read.table(paste0(loc_anndata,"/var.csv"),sep = ",",header = T)
obs<-read.table(paste0(loc_anndata,"/obs.csv"),sep = ",",header = T)
str(var)
str(obs)
colnames(X)<-var$index
rownames(X)<-obs$index

DimPlot(X)

resultTable<-read.table("GeneNames_EnsemblID.txt")
nrow(resultTable)


names_df<-resultTable[match(colnames(X),resultTable$wikigene_name),]
xXx<-X
colnames(xXx)<-names_df$ensembl_gene_id
dim(xXx)
xXx<-xXx[,c(which(is.na(colnames(xXx))==F))]
dim(xXx)
t(xXx)
write.table("Count_mtx_for_ChellphoneDB")
