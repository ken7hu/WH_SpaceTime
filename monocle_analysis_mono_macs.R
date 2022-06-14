#code for monocle analysis of monomac subsets
#for the mhcii lo cells
#redo the dimensionreduction for the mhcii lo subset of mono/mac cells
mm.sub.1 = SubsetData(mm_clean,ident.use = c('Mono_1','Mono_2','Mono_Mac_1','Mono_Mac_2','Mono_Mac_3','Mono_Mac_4','Mono_Mac_5','Mono_Mac_6'))
DefaultAssay(mm.sub.1)<-'SCT'
mm.sub.1 <- RunPCA(mm.sub.1,verbose = FALSE)
ElbowPlot(mm.sub.1,reduction = "pca",ndims = 50)
mm.sub.1 <- mm.sub.1 %>% RunUMAP( dims = 1:17, verbose = T,min.dist  = 0.03,n.neighbors = 40,seed.use = 12345)

DefaultAssay(mm.sub.1)<-'RNA'
cds <- as.cell_data_set(mm.sub.1)
cds <- cluster_cells(cds)
plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
integrated.sub <- subset(as.Seurat(cds), monocle3_partitions == 1)
cds <- as.cell_data_set(integrated.sub)
cds <- learn_graph(cds)
plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)
cds <- order_cells(cds)
plot_cells(
  cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE
)
#remove stray cluster of mito hi cells
macs_sub_part1 = subset(mm.sub.1,cells = colnames(integrated.sub))
DimPlot(macs_sub_part1)
cds <- as.cell_data_set(macs_sub_part1)
cds <- cluster_cells(cds)
plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
cds <- learn_graph(cds)
plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)
cds <- order_cells(cds)
p = plot_cells(
  cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE,
  cell_size = 0.7
)

cds <- estimate_size_factors(cds)
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(macs_sub_part1[["RNA"]])

#### plotting genes based on pseudotime
AFD_genes <- c("Pparg")

AFD_lineage_cds <- cds[rowData(cds)@rownames %in% AFD_genes,]

trajectories = plot_genes_in_pseudotime(AFD_lineage_cds,
                                        color_cells_by="Day", trend_formula = "~ splines::ns(pseudotime, df=7)",
                                        min_expr=1,label_by_short_name = FALSE,cell_size = 1.5,vertical_jitter = 0.01)
##### extract pseudotime stat for seurat usage
sudo = cds@principal_graph_aux@listData$UMAP$pseudotime
mm.sub.1<-AddMetaData(mm.sub.1,metadata = sudo,col.name = 'pseudotime')
v1 = VlnPlot(mm.sub.1,features = 'pseudotime',group.by = 'Day')
v2 = VlnPlot(mm.sub.1,features = 'pseudotime')
v2 = v2+scale_fill_manual(values = palette.use[1:8])
v1+v2+plot_layout(ncol=2)
#######FINAL PLOTTING##################
library(scales)
palette.use = c(hue_pal()(11)[c(1,3,5,7,9,11)],hue_pal()(11)[c(2,4,6,8,10)])
g = DimPlot(mm.sub.1,label = T)
g = g+theme(legend.position = 'None')+scale_color_manual(values = palette.use[1:8])
library(patchwork)
g +p+plot_layout(ncol =2)

##############################################
#Repeat for MHCII Hi Subset of Mono_Macs
mm.sub.2 = SubsetData(mm_clean,ident.use = c('Mono_MHCII','Mono_Mac_MHCII','Mono_Mac_MHCII_Mgl2'))
DefaultAssay(mm.sub.2)<-'SCT'
mm.sub.2 <- RunPCA(mm.sub.2,verbose = FALSE)
ElbowPlot(mm.sub.2,reduction = "pca",ndims = 50)
mm.sub.2 <- mm.sub.2 %>% RunUMAP( dims = 1:15, verbose = TRUE,min.dist = 0.15,n.neighbors = 40,seed.use = 12345,spread=0.5)
DefaultAssay(mm.sub.2)<-'RNA'
cds <- as.cell_data_set(mm.sub.2)
cds <- cluster_cells(cds)
plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
integrated.sub <- subset(as.Seurat(cds), monocle3_partitions == 1)
cds <- as.cell_data_set(integrated.sub)
cds <- learn_graph(cds)
plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)
cds <- order_cells(cds)
plot_cells(
  cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE
)
##### extract pseudotime stat for seurat usage
sudo = cds@principal_graph_aux@listData$UMAP$pseudotime
mm.sub.2<-AddMetaData(mm.sub.2,metadata = sudo,col.name = 'pseudotime')
v1 = VlnPlot(mm.sub.2,features = 'pseudotime',group.by = 'Day')
v2 = VlnPlot(mm.sub.1,features = 'pseudotime')
v2 = v2+scale_fill_manual(values = palette.use[1:8])
v1+v2+plot_layout(ncol=2)
#####subset out the stray clusters/;
cds.sub = cds[,clusters(cds, reduction_method = "UMAP")==1]
mm.sub.2.sub = mm.sub.2[,clusters(cds, reduction_method = "UMAP")==1]
#######
sudo = cds.sub@principal_graph_aux@listData$UMAP$pseudotime
mm.sub.2.sub<-AddMetaData(mm.sub.2.sub,metadata = sudo,col.name = 'pseudotime')
v1 = VlnPlot(mm.sub.2.sub,features = 'pseudotime',group.by = 'Day')
v2 = VlnPlot(mm.sub.2.sub,features = 'pseudotime')
v2 = v2+scale_fill_manual(values = palette.use[9:11])
v1+v2+plot_layout(ncol=2)
#####final plotting##############
palette.use = c(hue_pal()(11)[c(1,3,5,7,9,11)],hue_pal()(11)[c(2,4,6,8,10)])
g = DimPlot(mm.sub.2.sub,label = T)
g = g+theme(legend.position = 'None')+scale_color_manual(values = palette.use[9:11])
p = plot_cells(
  cds.sub,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE,
  cell_size = 0.7
)
g+p+plot_layout(ncol = 2)

