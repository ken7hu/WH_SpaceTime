library(readr)
library(stringr)
#read exported cytof csv
WH010_phemd <- read_csv("C:/Users/KHu/Downloads/WH013_phemd.csv")
#we had 40 channels
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

# loading into seurat -----------------------------------------------------
#load into seurat
pbmc10k.data <- Read10X(data.dir = "~/Krummel lab/phemd/WH013/gzipped/")
pbmc10k <- CreateSeuratObject(counts = pbmc10k.data)
# Set an Assay slot directly
#featch the raw values from the RNA assay
count.data <- GetAssayData(object = pbmc10k[["RNA"]], slot = "counts")
#arcsinh transform and create a new assay called 'CYTOF'
count.sinh = asinh((count.data@x))
count.data@x=count.sinh
pbmc10k[["CYTOF"]] <- CreateAssayObject(count.data)
pbmc10k <- SetAssayData(pbmc10k, assay = "CYTOF", slot = "scale.data", new.data = as.matrix(count.data))
all.genes = row.names(pbmc10k)
DefaultAssay(pbmc10k)<-'CYTOF'
pbmc10k <- RunPCA(pbmc10k, features = all.genes,assay = 'CYTOF')
#run clustering directly on the features
pbmc10k <- FindNeighbors(pbmc10k, features = all.genes,assay = 'CYTOF')
#pbmc_10k <- FindNeighbors(pbmc10k, dims = 1:40,reduction = "pca")
pbmc10k <- FindClusters(pbmc10k, resolution = 2)
pbmc10k <- RunUMAP(pbmc10k, features = all.genes,min.dist = 0.1,n.neighbors = 200,slot = "scale.data",verbose=TRUE,assay = 'CYTOF')
#assign the sample column from the cytof export
pbmc10k$sample = as.factor(WH010_phemd$sample)
levels(pbmc10k$sample)<-c('UW_1','UW_2','UW_3','D0.25_1','D0.25_2','D0.25_3','D01_1','D01_2','D01_3','D03_1','D03_2','D03_3',
                          'D05_1','D05_2','D05_3','D07_1','D07_2','D07_3','D10_1','D10_2','D10_3',
                          'D14_1','D14_2','D14_3','D21_1','D21_2','D21_3')
#for plotting
DimPlot(pbmc10k,label=TRUE)
# convert to phemd object -------------------------------------------------
# make a list of expression data for each sample: 
expression_matrix <- GetAssayData(pbmc10k, slot='data',assay = 'CYTOF')
expression_matrix <- as.matrix(expression_matrix)
expression_list <- list()
for(sample in unique(pbmc10k$sample)){
  print(sample)
  expression_list[[sample]] <- as.matrix(t(expression_matrix[,pbmc10k$sample == sample]))
}
phemd_obj <- createDataObj(expression_list, all.genes, names(expression_list))
phemd_obj <- bindSeuratObj(phemd_obj, pbmc10k)
batchIDs(phemd_obj) <- names(expression_list)
phemd_obj@seurat_obj@assays$RNA=phemd_obj@seurat_obj@assays$CYTOF
phemd_obj@seurat_obj$plt = phemd_obj@seurat_obj$sample
phemd_obj <- clusterIndividualSamples(phemd_obj, cell_model='seurat')
phemd_obj <- generateGDM(phemd_obj, cell_model='seurat', expn_type='umap', ndim=2)
emd_distmat <- compareSamples(phemd_obj)
group_assignments <- groupSamples(emd_distmat, distfun='hclust', ncluster=5)
phemd_dmap <- plotGroupedSamplesDmap(emd_distmat, group_assignments, pt_sz=1.5,n_dim = 2)
eigvects = phemd_dmap@eigenvectors
diff.df = as.data.frame(eigvects[,c(1,2)])

####for 2D plotting of WH010
days_list<-c('UW','UW','UW','D0.25','D0.25','D0.25','D01','D01','D01','D03','D03','D03',
                          'D05','D05','D05','D07','D07','D07','D10','D10','D10',
                          'D14','D14','D14','D21','D21','D21')
diff.df$Day = factor(days_list)
ggplot(data = diff.df,aes(x=DC1,y=-DC2,color = Day))+geom_point(size=2)+ theme_classic()+labs(color = 'Days Post-Wound')+
  theme(axis.title = element_text(size=20),axis.text = element_blank(),legend.text = element_text(size=16),legend.title = element_text(size=16))

#####to plot histograms of the various clusters by sanmple
diff.df$Day = sub("\\_.*", "", diff.df$Punch_ID)

pretty= c("red4","red2","darkorange1","goldenrod2","yellowgreen","chartreuse","darkgreen","skyblue","navy","darkorchid3","gray57")
sample.cellfreqs <- getSampleHistsByCluster(phemd_obj, group_assignments, cell_model='seurat')
plotSummaryHistograms(phemd_obj, group_assignments, cell_model='seurat', 
                      ncol.plot = 5, ax.lab.sz=1.3, title.sz=2)
