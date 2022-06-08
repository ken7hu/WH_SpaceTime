library(readr)
library(stringr)
WH013_phemd <- read_csv("C:/Users/KHu/Downloads/WH013_phemd.csv")
test = WH013_phemd[,1:40]
colnames(test)[1:40]<-str_extract(colnames(test)[1:40],"(?<=_).*")
features = data.frame(colnames(test)[1:40])
features$gene = (colnames(test)[1:40])
features$type = rep("Gene Expression",nrow(features))
#write mock features and barcodes .tsv files
write.table(features,col.names = FALSE,row.names = FALSE,sep="\t",
            file = "~/Krummel lab/phemd/WH013/features.tsv",quote = FALSE)
barcodes = c(1:nrow(test))
write.table(barcodes,col.names = FALSE,row.names = FALSE,sep="\t",
            file = "~/Krummel lab/phemd/WH013/barcodes.tsv",quote = FALSE)
library(readr)
#write a mtx sparse matrix
m<-t(as.matrix(test[,1:40]))
colnames(m)<-NULL
m <- Matrix(m, sparse = TRUE)
m <- as(m, "dgTMatrix")
writeMM(m,"~/Krummel lab/phemd/WH013/matrix.mtx")
#load into seurat
pbmc10k.data <- Read10X(data.dir = "~/Krummel lab/phemd/WH013/gzipped/")

pbmc10k <- CreateSeuratObject(counts = pbmc10k.data)

# Set an Assay slot directly
#featch the raw values from the RNA assay
count.data <- GetAssayData(object = pbmc10k[["RNA"]], slot = "counts")
#arcsinh transform and create a new assay called 'CYTOF'
blah2 = asinh((count.data@x))
count.data@x=blah2
pbmc10k[["CYTOF"]] <- CreateAssayObject(count.data)
pbmc10k <- SetAssayData(pbmc10k, assay = "CYTOF", slot = "scale.data", new.data = as.matrix(count.data))

all.genes = row.names(pbmc10k)
DefaultAssay(pbmc10k)<-'CYTOF'
pbmc_10k <- RunPCA(pbmc_10k, features = all.genes,assay = 'CYTOF')
#run clustering directly on the features
pbmc_10k <- FindNeighbors(pbmc_10k, features = all.genes,assay = 'CYTOF')
#pbmc_10k <- FindNeighbors(pbmc10k, dims = 1:40,reduction = "pca")
pbmc_10k <- FindClusters(pbmc_10k, resolution = 1)
pbmc_10k <- RunUMAP(pbmc_10k, features = all.genes,min.dist = 0.1,n.neighbors = 200,slot = "scale.data",verbose=TRUE,assay = 'CYTOF')
pbmc10k <- RunUMAP(pbmc10k, dims = 1:40,min.dist = 0.1,n.neighbors = 200,slot = "scale.data")
#for plotting
DimPlot(pbmc10k,label=TRUE)
DimPlot(pbmc10k,group.by = "sample")
FeaturePlot(pbmc10k,features = c("CD3e","CD4","CD8"))
#####################################################################
#move on to actual PhEMD
#######################################################################
#grab the asinh transformed cytof assay 
expression_matrix <- GetAssayData(pbmc10k, slot='counts',assay = 'CYTOF')
expression_matrix <- as.matrix(expression_matrix)
expression_list <- list()
pbmc10k$sample = pbmc10k$ADT_maxID
for(sample in unique(pbmc10k$sample)){
  print(sample)
  expression_list[[sample]] <- as.matrix(t(expression_matrix[,pbmc10k$sample == sample]))
}

#create a phemd object
phemd_obj <- createDataObj(expression_list,rownames(pbmc10k@assays$RNA) , names(expression_list))
phemd_obj <- bindSeuratObj(phemd_obj, pbmc21)
batchIDs(phemd_obj) <- names(expression_list)

phemd_obj@seurat_obj$plt = phemd_obj@seurat_obj$sample
phemd_obj <- clusterIndividualSamples(phemd_obj, cell_model='seurat')
phemd_obj <- generateGDM(phemd_obj, cell_model='seurat', expn_type='umap', ndim=2)
#calculate the earh movers distances
emd_distmat <- compareSamples(phemd_obj)
group_assignments <- groupSamples(emd_distmat, distfun='hclust', ncluster=2)
phemd_dmap <- plotGroupedSamplesDmap(emd_distmat, group_assignments, pt_sz=1.5)
phemd_dmap <- plotGroupedSamplesDmap(emd_distmat, cond, pt_sz=1.5)
#color scheme selsection
pretty= c("red4","red2","darkorange1","goldenrod2","yellowgreen","chartreuse","darkgreen","skyblue","navy","darkorchid3","gray57")
emd_distmat <- compareSamples(phemd_obj)
group_assignments <- groupSamples(emd_distmat, distfun='hclust', ncluster=5)
phemd_dmap <- plotGroupedSamplesDmap(emd_distmat, group_assignments, pt_sz=1.5)
phemd_dmap <- plotGroupedSamplesDmap(emd_distmat, cond, pt_sz=1.5)
phemd_dmap <- plotGroupedSamplesDmap(emd_distmat, group_assignments, pt_sz=1.5,n_dim = 2)
phemd_dmap <- plotGroupedSamplesDmap(emd_distmat, cond, pt_sz=1.5,n_dim=2)
#for 2D phemd
phemd_dmap <- plotGroupedSamplesDmap(emd_distmat, day_list, pt_sz=1.5,pt_label = names(expression_list),n_dim = 2)
eigvects = phemd_dmap@eigenvectors
diff.df = as.data.frame(eigvects[,c(1,2)])
diff.df$cond = cond
diff.df$Day = as.factor(day_list)
diff.df$Punch_ID = names(expression_list)
ggplot(data = diff.df,aes(x=DC1,y=DC2,color = Day))+geom_point(size=2)+
  geom_text_repel(aes(label=Day),color = "black",point.padding = 0.5,box.padding = 0.35,size=5)
ggplot(data = diff.df,aes(x=DC1,y=DC2,color = Day))+geom_point()+
  geom_text(aes(label=Day),color = "black")

####for 2D plotting of WH010 (figure 1E)
days_list = c('UW','D0.25','D1','D3','D5','D7','D10','D14','D21')
diff.df$Day = days_list[diff.df$cond]
diff.df$Day = factor(diff.df$Day,levels = days_list)


ggplot(data = diff.df,aes(x=DC1,y=DC2,color = Day))+geom_point(size=2)+ theme_classic()+labs(color = 'Days Post-Wound')+
  theme(axis.title = element_text(size=20),axis.text = element_blank(),legend.text = element_text(size=16),legend.title = element_text(size=16))

ggplot(data = diff.df,aes(x=DC1,y=DC2,color = Day))+geom_point(size=2)+ theme_classic()+labs(color = 'Days Post-Wound')+ 
  scale_color_viridis(discrete=TRUE,option = 'turbo') +
  theme(axis.title = element_text(size=20),axis.text = element_blank(),legend.text = element_text(size=16),legend.title = element_text(size=16))

#####for 2D plotting of WH012 (repeat of figure 1E)
diff.df$Day = sub("\\_.*", "", diff.df$Punch_ID)

pretty= c("red4","red2","darkorange1","goldenrod2","yellowgreen","chartreuse","darkgreen","skyblue","navy","darkorchid3","gray57")

