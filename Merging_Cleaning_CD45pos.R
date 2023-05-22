WH011_1<- NormalizeData(WH011_1)
WH011_2<- NormalizeData(WH011_2)
WH011_3<- NormalizeData(WH011_3)
WH011_4<- NormalizeData(WH011_4)
WH011_5<- NormalizeData(WH011_5)
pbmc.big <- merge(WH011_1, y = c(WH011_2,WH011_3,WH011_4,WH011_5), add.cell.ids = c("1","2","3","4","5"), project = "WH011",merge.data = TRUE)
pbmc.big$library = as.factor(as.integer((substr(colnames(pbmc.big),1,1))))
DefaultAssay(pbmc.big)<-'RNA'
pbmc.big <- PercentageFeatureSet(pbmc.big, pattern = "^mt-", col.name = "percent.mt")
pbmc.big <-subset(pbmc.big,subset = percent.mt<20)
pbmc.big<- SCTransform(pbmc.big, verbose = TRUE)
pbmc.big <- RunPCA(pbmc.big,verbose = FALSE)
ElbowPlot(pbmc.big,reduction = "pca")
pbmc.big <- FindNeighbors(pbmc.big, dims = 1:19, verbose = FALSE)
pbmc.big <- FindClusters(pbmc.big, verbose = FALSE,resolution = 0.8)
pbmc.big <- pbmc.big %>% RunUMAP( dims = 1:19, verbose = TRUE)
DimPlot(pbmc.big,label=TRUE)
FeaturePlot(pbmc.big,features = c('Ptprc','Col1a1'))
##remove clusters 8 and 6 for being too high in mito
pbmc.big<-subset(pbmc.big,idents = c(8,26),invert=TRUE)
pbmc.big <- RunPCA(pbmc.big,verbose = FALSE)
ElbowPlot(pbmc.big,reduction = "pca")
pbmc.big <- FindNeighbors(pbmc.big, dims = 1:19, verbose = FALSE)
pbmc.big <- FindClusters(pbmc.big, verbose = FALSE,resolution = 0.4)
pbmc.big <- pbmc.big %>% RunUMAP( dims = 1:19, verbose = TRUE)

DimPlot(pbmc.big,cells.highlight = colnames(pbmc.big)[pbmc.big$nCount_RNA<6250])
ggplot(blah, aes(fill=library, y=Freq, x=cluster)) + 
  geom_bar(position="fill", stat="identity")


#these clusters are enriched for libraries 4or 5
pbmc.big <- FindClusters(pbmc.big, verbose = FALSE,resolution = 1)
pbmc.big<-subset(pbmc.big,idents = c(1,2,3,5,7,10),invert=TRUE)
pbmc.sub<-subset(pbmc.big,subset = percent.mt<25)
#perform SC Transform
pbmc.sub<- SCTransform(pbmc.sub, vars.to.regress = "percent.mt", verbose = FALSE)
pbmc.sub <- RunPCA(pbmc.sub,verbose = FALSE)
ElbowPlot(pbmc.sub,reduction = "pca",ndims = 50)
pbmc.sub <- FindNeighbors(pbmc.sub, dims = 1:28, verbose = FALSE)
pbmc.sub <- FindClusters(pbmc.sub, verbose = FALSE,resolution = 0.8)
pbmc.sub <- pbmc.sub %>% RunUMAP( dims = 1:28, verbose = TRUE)
DimPlot(pbmc.sub,label=TRUE)
DimPlot(pbmc.sub,group.by = 'library')
blah = interaction(pbmc.big$library,Idents(pbmc.big)) %>% table() %>% as.data.frame()
blah$library = as.factor(rep(c(1:5),30))
blah$cluster = as.factor(floor(c(0:149)/5))
ggplot(blah, aes(fill=library, y=Freq, x=cluster)) + 
  geom_bar(position="fill", stat="identity")

pbmc.sub <- FindNeighbors(pbmc.sub, dims = 1:28, verbose = FALSE)
pbmc.sub <- FindClusters(pbmc.sub, verbose = FALSE,resolution = 0.6)
DefaultAssay(pbmc.sub)<-'RNA'
pbmc.markers <- FindAllMarkers(pbmc.sub, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.2)
b_plasma_harmony_0.2_Markers_padj0.1_Top5 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DotPlot(pbmc.sub,
        features = unique(b_plasma_harmony_0.2_Markers_padj0.1_Top5$gene),
        cols = "RdBu", assay = "RNA") + theme(axis.text.x=element_text(angle=45, hjust = 1, size = 10), axis.text.y=element_text(size = 10), text = element_text(size = 14)) + coord_flip()

###this is where it is saved as a rds under "20210309_wh_clean.rds")

# Final clean step --------------------------------------------------------

pbmc.sub2<-subset(pbmc.sub,subset = nFeature_RNA>250)
pbmc.sub2<-subset(pbmc.sub2,subset = nCount_RNA>400)
pbmc.sub2 <- RunPCA(pbmc.sub2,verbose = FALSE)
ElbowPlot(pbmc.sub2,reduction = "pca",ndims = 50)
pbmc.sub2 <- FindNeighbors(pbmc.sub2, dims = 1:28, verbose = FALSE)
pbmc.sub2 <- FindClusters(pbmc.sub2, verbose = FALSE,resolution = 0.5)
pbmc.sub2 <- pbmc.sub2 %>% RunUMAP( dims = 1:28, verbose = TRUE)
DimPlot(pbmc.sub2,label=TRUE)
#####split into cd45+/- pops
DefaultAssay(pbmc.sub2)<-'RNA'
pbmc.markers <- FindAllMarkers(pbmc.sub2, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.2)
b_plasma_harmony_0.2_Markers_padj0.1_Top5 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DotPlot(pbmc.sub2,
        features = unique(b_plasma_harmony_0.2_Markers_padj0.1_Top5$gene),
        cols = "RdBu", assay = "RNA") + theme(axis.text.x=element_text(angle=45, hjust = 1, size = 10), axis.text.y=element_text(size = 10), text = element_text(size = 14)) + coord_flip()
#clusters c(0,1,2,3,5,7,9,18,20,21) are all cd45+
cd45.pos = subset(pbmc.sub2,idents = c(0,1,2,3,5,7,9,18,20,21))

# CD45 positive drill down ------------------------------------------------

cd45.pos <- CellCycleScoring(cd45.pos, s.features = mouse.s.genes, g2m.features = mouse.g2m.genes, set.ident = FALSE)
cd45.pos<- SCTransform(cd45.pos, vars.to.regress = c("percent.mt"), verbose = TRUE)
cd45.pos <- RunPCA(cd45.pos,verbose = FALSE)
ElbowPlot(cd45.pos,reduction = "pca",ndims = 50)
DimHeatmap(cd45.pos,dims = c(9),slot = 'data')
cd45.pos <- FindNeighbors(cd45.pos, dims = 1:23, verbose = FALSE)
cd45.pos <- FindClusters(cd45.pos, verbose = FALSE,resolution = 0.9)
cd45.pos <- cd45.pos %>% RunUMAP( dims = 1:23, verbose = TRUE,min.dist = 0.1)
DimPlot(cd45.pos,label=TRUE)
table(Idents(cd45.pos))
DefaultAssay(cd45.pos)<-'RNA'
pbmc.markers <- FindAllMarkers(cd45.pos, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.2)
b_plasma_harmony_0.2_Markers_padj0.1_Top5 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DotPlot(cd45.pos,
        features = unique(b_plasma_harmony_0.2_Markers_padj0.1_Top5$gene),
        cols = "RdBu", assay = "RNA") + theme(axis.text.x=element_text(angle=45, hjust = 1, size = 10), axis.text.y=element_text(size = 10), text = element_text(size = 14)) + coord_flip()

ggplot(blah, aes(fill=type, y=Freq, x=Punch_ID)) + 
  geom_bar(position="fill", stat="identity",color="white")+
  theme(axis.text.x = element_text(angle = 45))

new.cluster.ids <- c("Mono_Mac_1", "Mono_Mac_mito_1", "Mono_Mac_2", "Mono", "Mono_Mac_MHCII", "Neut_1", 
                     "Collagen", "Neuts_2", "Mono_Mac_4","Neuts_3",
                     "Neuts_4","Mono_Mac_5","Mono_Mac_mito_2","T_1","Mono_Mac_6",
                     "DC","Mono_Vcan","TNK","Mono_Mac_Mgl2","T_2","Mast Cells",
                     "Prolif","B")
names(new.cluster.ids) <- levels(cd45.pos)
cd45.pos$stashIdent_0_9 = Idents(cd45.pos)
cd45.pos <- RenameIdents(cd45.pos, new.cluster.ids)

cluster.orders <- c("Mono_Mac_1",  "Mono_Mac_2", "Mono", "Mono_Mac_MHCII","Mono_Mac_4", "Mono_Mac_5","Mono_Mac_mito_1", 
                    "Mono_Mac_Mgl2","Mono_Mac_mito_2","Mono_Vcan","Mono_Mac_6","DC","Collagen", 
                    "Neut_1","Neuts_2", "Neuts_3", "Neuts_4",
                    "TNK","T_1","T_2","Mast Cells",
                    "Prolif","B")

blah_sub$type = factor(blah_sub$type,levels = cluster.orders)

######for replotting with monomacmito1 and collagen removed

cd45.pos = subset(cd45.pos,idents = c('Mono_MAC_mito_1','Collagen'),invert=TRUE)
cd45.pos <- cd45.pos %>% RunUMAP( dims = 1:23, verbose = TRUE,min.dist = 0.2)
DimPlot(cd45.pos,label = TRUE)
blah_sub$type = sub(".*\\.","",blah_sub$.)
blah_sub$Punch_ID = sub("\\..*", "", blah_sub$.)

###now CC regression
cd45.pos$annotation_1<-Idents(cd45.pos)
cd45.pos<- SCTransform(cd45.pos, vars.to.regress = c("percent.mt","S.Score","G2M.Score"), verbose = TRUE)
cd45.pos <- RunPCA(cd45.pos,verbose = FALSE)
ElbowPlot(cd45.pos,reduction = "pca",ndims = 50)
cd45.pos <- FindNeighbors(cd45.pos, dims = 1:28, verbose = FALSE)
cd45.pos <- FindClusters(cd45.pos, verbose = FALSE,resolution = 0.9)
cd45.pos <- cd45.pos %>% RunUMAP( dims = 1:28, verbose = TRUE,min.dist = 0.1)
DimPlot(cd45.pos,label = TRUE)

##end here for mac/T/neut subestting
#subset out mono/macs:
macs = subset(cd45.pos,idents = c(0,1,3,4,6,8,9,11,14))
neuts = subset(cd45.pos,idents = c(2,5,7))
t = subset(cd45.pos,idents = c(10,15,18))

###recluster and plot macs:
DefaultAssay(macs)<-'SCT'
macs<- SCTransform(macs, vars.to.regress = c("percent.mt","S.Score","G2M.Score"), verbose = TRUE)
macs <- RunPCA(macs,verbose = FALSE)
ElbowPlot(macs,reduction = "pca",ndims = 50)
macs <- FindNeighbors(macs, dims = 1:30, verbose = FALSE)
macs <- FindClusters(macs, verbose = FALSE,resolution = 0.6)
macs <- macs %>% RunUMAP( dims = 1:30, verbose = TRUE,min.dist = 0.1)
DimPlot(macs,label = TRUE)
DefaultAssay(macs)<-'RNA'
pbmc.markers <- FindAllMarkers(macs, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.2)
b_plasma_harmony_0.2_Markers_padj0.1_Top5 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DotPlot(macs,
        features = unique(b_plasma_harmony_0.2_Markers_padj0.1_Top5$gene),
        cols = "RdBu", assay = "RNA") + theme(axis.text.x=element_text(angle=45, hjust = 1, size = 10), axis.text.y=element_text(size = 10), text = element_text(size = 14)) + coord_flip()
#remove cluster 7, neut contamination? 12 is langerhans, 11 is more DC's, 14 and 15 other remaining doublets
macs<-subset(macs,idents = c(7,8,11,14,15),invert=TRUE)
macs <- RunPCA(macs,verbose = FALSE)
ElbowPlot(macs,ndims = 30)
DefaultAssay(macs)<-'SCT'
macs <- FindNeighbors(macs, dims = 1:12, verbose = FALSE)
macs <- FindClusters(macs, verbose = FALSE,resolution = 0.7)
macs <- macs %>% RunUMAP( dims = 1:12, verbose = TRUE,min.dist = 0.005)
DimPlot(macs,label = T)


