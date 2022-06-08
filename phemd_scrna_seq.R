#PHEMD code for going directly from Seurat scRNA-Seq experiment (i.e. Fig 2J)
# make a list of expression data for each sample: 
expression_matrix <- GetAssayData(cd45.pos, slot='data',assay = 'RNA')
expression_matrix <- as.matrix(expression_matrix)
expression_list <- list()
cd45.pos$sample = cd45.pos$Punch_ID
cd45.pos$annotated = Idents(cd45.pos)
Idents(cd45.pos)<-as.numeric(Idents(cd45.pos))
for(sample in unique(cd45.pos$sample)){
  print(sample)
  expression_list[[sample]] <- as.matrix(t(expression_matrix[,cd45.pos$sample == sample]))
}

phemd_obj <- createDataObj(expression_list,rownames(cd45.pos@assays$RNA) , names(expression_list))
phemd_obj <- bindSeuratObj(phemd_obj, cd45.pos)
batchIDs(phemd_obj) <- names(expression_list)
phemd_obj@seurat_obj$plt = phemd_obj@seurat_obj$sample
phemd_obj <- clusterIndividualSamples(phemd_obj, cell_model='seurat')
phemd_obj <- generateGDM(phemd_obj, cell_model='seurat', expn_type='umap', ndim=2)

emd_distmat <- compareSamples(phemd_obj)
group_assignments <- groupSamples(emd_distmat, distfun='hclust', ncluster=2)
phemd_dmap <- plotGroupedSamplesDmap(emd_distmat, group_assignments, pt_sz=1.5,n_dim = 2)
phemd_dmap <- plotGroupedSamplesDmap(emd_distmat, cond, pt_sz=1.5)
phemd_dmap <- plotGroupedSamplesDmap(emd_distmat, day_list, pt_sz=1.5,pt_label = names(expression_list),n_dim = 2)
pretty= c("red4","red2","darkorange1","goldenrod2","yellowgreen","chartreuse","darkgreen","skyblue","navy","darkorchid3","gray57")
eigvects = phemd_dmap@eigenvectors
diff.df = as.data.frame(eigvects[,c(1,2)])
diff.df$Day = as.factor(day_list)
diff.df$Punch_ID = names(expression_list)
ggplot(data = diff.df,aes(x=DC1,y=DC2,color = Day))+geom_point()+
  geom_text_repel(aes(label=Punch_ID),color = "black",point.padding = 0.5,box.padding = 0.35)

ggplot(data = diff.df,aes(x=DC1,y=DC2,color = Day))+geom_point()+
  geom_text(aes(label=Punch_ID),color = "black")

ggplot(data = diff.df,aes(x=DC1,y=DC2,color = Day))+geom_point(size=2)+theme_classic()+labs(color = 'Days Post-Wound')+ 
  theme(axis.title = element_text(size=20),axis.text = element_blank(),legend.text = element_text(size=16),legend.title = element_text(size=16))+
  geom_text_repel(aes(label=Punch_ID),color = "black",point.padding = 0.1,box.padding = 0.35,label.padding = 0.3,direction = "both",force = 2,size=5)


g1 = ggplot(data = diff.df,aes(x=DC1,y=DC2,color = Day))+geom_point(size=4)+theme_classic()+labs(color = 'Days Post-Wound')+ 
  theme(axis.title = element_text(size=20),axis.text = element_blank(),legend.text = element_text(size=16),legend.title = element_text(size=16))+
  geom_text_repel(aes(label=Punch_ID),color = "black",box.padding = 0.5,point.padding = 0.1,direction = "both",force = 5,size=7,segment.size = 0.75)
g1+geom_mark_hull(aes(fill = Day),concavity = 1,expand = unit(1.4, "mm"))
g1+geom_mark_hull(aes(fill = Day),concavity = 1,expand = unit(1.4, "mm"),radius = unit(1.4, "mm"))
g1+geom_mark_ellipse(aes(fill = Day))

###for umap dim red on the distance matrix emd_distmat
library(umap)
test_umap = umap(emd_distmat,input = "dist")
day_list = (c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7))
##Plot Summary Histograms
plotSummaryHistograms(phemd_obj, day_list, cell_model='seurat',  
                      ncol.plot = 5, ax.lab.sz=1.3, title.sz=2)
