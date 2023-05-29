###### integrated CD45neg object processing to identify cell subsets (Figure 3A and Supplemental Figures 3A + 3B) and subsetting out fibroblasts for further analysis (see separate Rscript file)
library(Seurat)
library(sctransform)
library(dplyr)
library(ggplot2)

#set working directory and load CD45-negative Seurat object
setwd("/Users/Downloads")
CD45neg <- readRDS("CD45neg.rds")

#add cell cycle scores to object; load cell cycle genes as csv file first (copied from human gene list)
library(readr)
mouse_g2m_genes <- read_csv("mouse_g2m_genes.csv")
#turn into vector type 'character'
mouse_g2m_genes <- mouse_g2m_genes$x
mouse_s_genes <- read_csv("mouse_s_genes.csv")
mouse_s_genes <- mouse_s_genes$x

CD45neg <- CellCycleScoring(CD45neg, s.features = mouse_s_genes, g2m.features = mouse_g2m_genes, set.ident = FALSE)
#SCT transform the object
CD45neg <- SCTransform(CD45neg, vars.to.regress = c("percent.mt", 'S.Score', 'G2M.Score'), verbose = T)

#run PCA and generate UMAP for plotting
CD45neg <- RunPCA(CD45neg, verbose = FALSE)
DimHeatmap(CD45neg, dims = 20:34, cells = 500, balanced = TRUE, nfeatures = 20)
ElbowPlot(CD45neg, ndims = 50)
CD45neg <- RunUMAP(CD45neg, dims = 1:21, n.neighbors = 100, min.dist = 0.3, spread = 1, a = NULL, b = NULL, verbose = FALSE, seed.use = 21212) #3 min
CD45neg <- FindNeighbors(CD45neg, dims = 1:21, k.param = 100, verbose = FALSE)
#FindClusters and plot UMAP
CD45neg_1.2 <- FindClusters(CD45neg, verbose = TRUE, algorithm = 1, resolution = 1.2, random.seed = 21212)  
DimPlot(CD45neg_1.2, reduction='umap', label = T) + labs(color = "1.2")

#FeaturePlotting to roughly identify subsets based on selected markers
DefaultAssay(CD45neg_1.2) <- 'RNA'
FeaturePlot(CD45neg_1.2, features = c('Ptprc', 'S100a8', 'H2-Ab1', 'Cd3e')) #immune
FeaturePlot(CD45neg_1.2, features = c('Col1a1', 'Pecam1', 'Acta2', 'Krt5')) #Fibro, endo, myofibro, epi
FeaturePlot(CD45neg_1.2, features = c('Fgf7', 'Hhip', 'Mlana', 'Kit')) #hair, melano
FeaturePlot(CD45neg_1.2, features = c('Krt17', 'Krt79', 'Cd34', 'Lgr5')) #hair follicle, HFSC
FeaturePlot(CD45neg_1.2, features = c('Krt5', 'Lgals7', 'Ly6d')) #epidermal

#unbiased marker finding
CD45neg_1.2_Markers <- FindAllMarkers(CD45neg_1.2,
                                           test.use='poisson',
                                           only.pos=TRUE,
                                           min.pct=0.25,
                                           logfc.threshold=0.25,
                                           assay='RNA')
CD45neg_1.2_Markers_padj0.1 <- CD45neg_1.2_Markers[which(CD45neg_1.2_Markers$p_val_adj<0.1),]
CD45neg_1.2_Markers_padj0.1 <- CD45neg_1.2_Markers_padj0.1[order(CD45neg_1.2_Markers_padj0.1$avg_log2FC,decreasing = TRUE),]
CD45neg_1.2_Markers_padj0.1 <- CD45neg_1.2_Markers_padj0.1[order(CD45neg_1.2_Markers_padj0.1$cluster,decreasing = FALSE),]
CD45neg_1.2_Markers_padj0.1_Top10 <- CD45neg_1.2_Markers_padj0.1 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
pdf("CD45neg_1.2_BubblePlot_RNA.pdf", width = 7, height = 7)
DotPlot(CD45neg_1.2,
        features = unique(CD45neg_1.2_Markers_padj0.1_Top10$gene),
        cols = 'RdBu', assay = "RNA") + theme(axis.text.x=element_text(angle=45, hjust = 1, size = 10), axis.text.y=element_text(size = 4), text = element_text(size = 14)) + coord_flip()
dev.off()

#barchart of all culsters by sequencing library
y <- table(CD45neg_1.2@meta.data[,c('library', 'seurat_clusters')])
y <- data.frame(y)
p <- ggplot(data=y, aes(fill = library, x=seurat_clusters, y=Freq))
p + geom_bar(position='fill', stat='identity')

#subset CD45neg object to remove clusteres 0, 1, 2, 3 (they predominantly come from library 4 & 5 and have myeloid marker genes) -> clean out potential doublets
CD45neg <- subset(CD45neg_1.2, idents = c('0', '1', '2', '3'), invert = TRUE)

#re-analyze subsetted, clean object
CD45neg <- CD45neg %>% SCTransform(vars.to.regress = c("percent.mt", 'S.Score', 'G2M.Score'), verbose = FALSE)
CD45neg <- RunPCA(CD45neg, verbose = FALSE)
DimHeatmap(CD45neg, dims = 20:34, cells = 500, balanced = TRUE, nfeatures = 20)
CD45neg <- RunUMAP(CD45neg, dims = 1:20, n.neighbors = 100, min.dist = 0.3, spread = 1, a = NULL, b = NULL, verbose = FALSE, seed.use = 21212)
CD45neg <- FindNeighbors(CD45neg, dims = 1:20, k.param = 100, verbose = FALSE)
#FindClusters and plot UMAP
CD45neg_2 <- FindClusters(CD45neg, verbose = TRUE, algorithm = 1, resolution = 2, random.seed = 212112)  
DimPlot(CD45neg_2, reduction='umap', label = T) + labs(color = "2")
#identify DEGs
CD45neg_2_Markers <- FindAllMarkers(CD45neg_2,
                                      test.use='poisson',
                                      only.pos=TRUE,
                                      min.pct=0.25,
                                      logfc.threshold=0.25,
                                      assay='RNA')
CD45neg_2_Markers_padj0.1 <- CD45neg_2_Markers[which(CD45neg_2_Markers$p_val_adj<0.1),]
CD45neg_2_Markers_padj0.1 <- CD45neg_2_Markers_padj0.1[order(CD45neg_2_Markers_padj0.1$avg_log2FC,decreasing = TRUE),]
CD45neg_2_Markers_padj0.1 <- CD45neg_2_Markers_padj0.1[order(CD45neg_2_Markers_padj0.1$cluster,decreasing = FALSE),]
CD45neg_2_Markers_padj0.1_Top5 <- CD45neg_2_Markers_padj0.1 %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
#plot DEGs in DotPlot format
DotPlot(CD45neg_2,
        features = unique(CD45neg_2_Markers_padj0.1_Top5$gene),
        cols = 'RdBu', assay = "RNA") + theme(axis.text.x=element_text(angle=45, hjust = 1, size = 10), axis.text.y=element_text(size = 8), text = element_text(size = 14))

#rename & reorder clusters based on DEGs
annotations <- c(
  "fibroblast_1", #0
  "keratinocyte_1", #1
  "fibroblast_5", #2
  "endothelial_cell_1", #3
  "endothelial_cell_2", #4
  "keratinocyte_2", #5
  "dermal_sheath_papilla", #6
  "endothelial_cell_3", #7 
  "fibroblast_2", #8 
  "keratinocyte_3", #9
  "fibroblast_3", #10
  "keratinocyte_4", #11
  "fibroblast_4", #12
  "immune", #13
  "melanocytes", #14
  "vascular_smooth_muscle_1", #15
  "sebaceous_gland", #16
  "keratinocyte_5", #17
  "vascular_smooth_muscle_2" #18
)
names(annotations) <- levels(CD45neg_2)
CD45neg_2 <- RenameIdents(CD45neg_2, annotations)
CD45neg_2@meta.data$type <- Idents(CD45neg_2)
DimPlot(CD45neg_2, reduction = 'umap', label = T) + NoLegend()
#reorder plots
levels(CD45neg_2) <- c("fibroblast_1", "fibroblast_2", "fibroblast_3", "fibroblast_4", "fibroblast_5", 
                       "keratinocyte_1", "keratinocyte_2", "keratinocyte_3", "keratinocyte_4", "keratinocyte_5", 
                       "endothelial_cell_1", "endothelial_cell_2", "endothelial_cell_3", "sebaceous_gland", "melanocytes", 
                       "vascular_smooth_muscle_1", "vascular_smooth_muscle_2", "dermal_sheath_papilla", "immune")

#plot Figure S3A
DotPlot(CD45neg_2,
        features = unique(CD45neg_2_Markers_padj0.1_Top5$gene),
        cols = 'RdBu', assay = "RNA") + theme(axis.text.x=element_text(angle=45, hjust = 1, size = 10), axis.text.y=element_text(size = 8), text = element_text(size = 14))

#make TilePlots for each cell subset -> Figure S3B
#prep data.frame:
#generate data.frame 'y' with one column of all clusters and one column of all PunchID
y <- table(CD45neg_2@meta.data[,c('type', 'PunchID')])
y <- data.frame(y)
#make temporary data.frame with # of each cell per PunchID
temp <- table(CD45neg_2@meta.data$PunchID)
temp <- data.frame(temp)
#add those numbers to y data.frame; matched by PunchID
y$total <- temp$Freq[match(y$PunchID, temp$Var1)]
#create column colled percents in y that is frequency of cell / PunchID
y$percents <- (y$Freq/y$total)*100
#add columns Day and Space with sub function; '.*' means keep everything before/after
y$Day = sub("_.*","",y$PunchID)
y$Space = sub(".*_","",y$PunchID)
#call data.frame blah
blah <- y

library(gridExtra)
plotlist = list()
for (i in 1:length(unique(blah$type))){
  df = subset(blah,subset = type == unique(blah$type)[i])
  uw_temp = df$percents[df$Day=='UW']
  df = subset(df,subset = Day != 'UW',invert=TRUE)
  df$percents_diff = df$percents-uw_temp
  rm(p1)
  p1 <-ggplot(df, aes(x = Day, y=Space)) +
    geom_tile(aes(fill = percents_diff),color='black') +
    geom_text(aes(label = round(percents, 1)))+
    scale_fill_gradient2(low = 'blue',mid = 'white',high = 'red',midpoint = 0)+
    labs(title = paste(unique(blah$variable)[i],'UW % :',as.character(round(uw_temp,3)),sep = ' '))+
    coord_fixed()+
    theme_void()+
    theme(axis.text = element_text(size = 16),legend.title = element_blank())
  plotlist[[i]]<-p1
}
##### TilePLot all at once
grid.arrange(grobs = plotlist, ncol = 5) 

##### TilePlot cluster individually
#exchange subset = type == 'cell_type' for cell_type of interest
df = subset(blah,subset = type == 'keratinocyte_1')
uw_temp = df$percents[df$Day=='UW']
df = subset(df,subset = Day != 'UW',invert=TRUE)
df$percents_diff = df$percents-uw_temp
g1 <-ggplot(df, aes(x = Day, y=Space)) +
  geom_tile(aes(fill = percents_diff),color='black') +
  geom_text(aes(label = round(percents, 1)), size=6)+
  scale_fill_gradient2(low = 'blue',mid = 'white',high = 'red',midpoint = 0)+
  labs(title = paste('keratinocyte_1','UW % :',as.character(round(uw_temp,1)),sep = ' '))+
  coord_fixed()+
  theme_void()+
  theme(axis.text = element_text(size = 16),legend.title = element_blank())
g1
df = subset(blah,subset = type == 'keratinocyte_2')
uw_temp = df$percents[df$Day=='UW']
df = subset(df,subset = Day != 'UW',invert=TRUE)
df$percents_diff = df$percents-uw_temp
g2 <-ggplot(df, aes(x = Day, y=Space)) +
  geom_tile(aes(fill = percents_diff),color='black') +
  geom_text(aes(label = round(percents, 1)), size=6)+
  scale_fill_gradient2(low = 'blue',mid = 'white',high = 'red',midpoint = 0)+
  labs(title = paste('keratinocyte_2','UW % :',as.character(round(uw_temp,1)),sep = ' '))+
  coord_fixed()+
  theme_void()+
  theme(axis.text = element_text(size = 16),legend.title = element_blank())
g2
df = subset(blah,subset = type == 'keratinocyte_3')
uw_temp = df$percents[df$Day=='UW']
df = subset(df,subset = Day != 'UW',invert=TRUE)
df$percents_diff = df$percents-uw_temp
g3 <-ggplot(df, aes(x = Day, y=Space)) +
  geom_tile(aes(fill = percents_diff),color='black') +
  geom_text(aes(label = round(percents, 1)), size=6)+
  scale_fill_gradient2(low = 'blue',mid = 'white',high = 'red',midpoint = 0)+
  labs(title = paste('keratinocyte_3','UW % :',as.character(round(uw_temp,1)),sep = ' '))+
  coord_fixed()+
  theme_void()+
  theme(axis.text = element_text(size = 16),legend.title = element_blank())
g3
df = subset(blah,subset = type == 'keratinocyte_4')
uw_temp = df$percents[df$Day=='UW']
df = subset(df,subset = Day != 'UW',invert=TRUE)
df$percents_diff = df$percents-uw_temp
g4 <-ggplot(df, aes(x = Day, y=Space)) +
  geom_tile(aes(fill = percents_diff),color='black') +
  geom_text(aes(label = round(percents, 1)), size=6)+
  scale_fill_gradient2(low = 'blue',mid = 'white',high = 'red',midpoint = 0)+
  labs(title = paste('keratinocyte_4','UW % :',as.character(round(uw_temp,1)),sep = ' '))+
  coord_fixed()+
  theme_void()+
  theme(axis.text = element_text(size = 16),legend.title = element_blank())
g4
df = subset(blah,subset = type == 'keratinocyte_5')
uw_temp = df$percents[df$Day=='UW']
df = subset(df,subset = Day != 'UW',invert=TRUE)
df$percents_diff = df$percents-uw_temp
g5 <-ggplot(df, aes(x = Day, y=Space)) +
  geom_tile(aes(fill = percents_diff),color='black') +
  geom_text(aes(label = round(percents, 1)), size=6)+
  scale_fill_gradient2(low = 'blue',mid = 'white',high = 'red',midpoint = 0)+
  labs(title = paste('keratinocyte_5','UW % :',as.character(round(uw_temp,1)),sep = ' '))+
  coord_fixed()+
  theme_void()+
  theme(axis.text = element_text(size = 16),legend.title = element_blank())
g5
(g1 | g2 | g3 | g4 | g5)
