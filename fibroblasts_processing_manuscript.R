library(Seurat)
library(sctransform)
library(dplyr)
library(ggplot2)

#set working directory and load CD45-negative Seurat object
setwd("/Users/nickkuhn/Downloads/spunchseq")
CD45neg_2 <-readRDS("CD45neg_2.rds") #read in processed CD45neg object

CD45neg_2$annotation_1 = Idents(CD45neg_2)
fibros <- subset(CD45neg_2, idents = c('fibroblast_1', 'fibroblast_2', 'fibroblast_3', 'fibroblast_4', 'fibroblast_5'))
fibros <- FindVariableFeatures(fibros,assay = 'SCT')
fibros <- RunPCA(fibros,verbose = FALSE)
ElbowPlot(fibros,reduction = "pca",ndims = 40)
fibros <- fibros %>% RunUMAP( dims = 1:12, verbose = FALSE,seed.use = 12345,n.neighbors = 30,spread=0.9)

#rename clusters
annotations <- c(
  "Fibro_1", #fibroblast_1
  "Fibro_2", #fibroblast_2
  "Fibro_3", #fibroblast_3
  "Fibro_4", #fibroblast_4
  "Fibro_5" #fibroblast_5
 )
names(annotations) <- levels(fibros)
fibros <- RenameIdents(fibros, annotations)
fibros@meta.data$type <- Idents(fibros)
#Plot UMAP -> Figure 3B
DimPlot(fibros, reduction = 'umap', label = T) + NoLegend()

#Run FindAllMarkers to identify DEGS
Idents(fibros) <- 'type'
fibro_Markers <- FindAllMarkers(fibros,
                                test.use='poisson',
                                only.pos=TRUE,
                                min.pct=0.25,
                                logfc.threshold=0.25,
                                assay='RNA')
fibro_Markers_padj0.1 <- fibro_Markers[which(fibro_Markers$p_val_adj<0.1),]
fibro_Markers_padj0.1 <- fibro_Markers_padj0.1[order(fibro_Markers_padj0.1$avg_log2FC,decreasing = TRUE),]
fibro_Markers_padj0.1 <- fibro_Markers_padj0.1[order(fibro_Markers_padj0.1$cluster,decreasing = FALSE),]

#plot DEG DotPlot -> Figure S3C
fibro_Markers_padj0.1_Top5 <- fibro_Markers_padj0.1 %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
pdf("fibros_DotPlot_RNA.pdf", width = 6, height = 5)
DotPlot(fibros,
        features = unique(fibro_Markers_padj0.1_Top5$gene),
        cols = 'RdBu', assay = "RNA") + theme(axis.text.x=element_text(angle=45, hjust = 1, size = 10), axis.text.y=element_text(size = 10), text = element_text(size = 14)) + coord_flip()
dev.off()

#plot VlnPlots for 'iCAF', myofibroblast, 'universal' fibroblast markers -> Figure S3D
s1 <- VlnPlot(fibros, c('Il1b', 'Cxcl2', 'Ccl2'), stack = T, fill.by = 'ident', flip = T) + NoLegend() + theme(axis.title.x = element_blank())
s2 <- VlnPlot(fibros, c('Acta2', 'Tagln', 'Mmp11'), stack = T, fill.by = 'ident', flip = T) + NoLegend() + theme(axis.title.x = element_blank())
s3 <- VlnPlot(fibros, c('Col15a1', 'Pi16', 'Cd34'), stack = T, fill.by = 'ident', flip = T) + NoLegend() + theme(axis.title.x = element_blank())
pdf('fb_stacked_vlnplots.pdf', width = 8, height = 3.5)
s1 | s2 | s3
dev.off()

#DimPlot highlight cells by day -> Figure 3C
gg1 <- DimPlot(fibros, pt.size = 0.1, cells.highlight = colnames(fibros)[fibros$Day=='D00'], sizes.highlight = 0.1) +theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) + NoLegend()
gg2 <- DimPlot(fibros, pt.size = 0.1, cells.highlight = colnames(fibros)[fibros$Day=='D01'], sizes.highlight = 0.1) +theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) + NoLegend()
gg3 <- DimPlot(fibros, pt.size = 0.1, cells.highlight = colnames(fibros)[fibros$Day=='D03'], sizes.highlight = 0.1) +theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) + NoLegend()
gg4 <- DimPlot(fibros, pt.size = 0.1, cells.highlight = colnames(fibros)[fibros$Day=='D07'], sizes.highlight = 0.1) +theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) + NoLegend()
gg5 <- DimPlot(fibros, pt.size = 0.1, cells.highlight = colnames(fibros)[fibros$Day=='D14'], sizes.highlight = 0.1) +theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) + NoLegend()
pdf('fibro_dimplot_by_day2.pdf', width = 2, height = 8.6)
(gg1 / gg2 / gg3 / gg4 / gg5)
dev.off()

#plot cell cluster % / day -> Figure 3D (left)
Idents(fibros) <- 'type'
y <- table(fibros@meta.data[,c('type', 'Day')])
y <- data.frame(y)
#make temporary data.frame with # of each cell per Day
temp <- table(fibros@meta.data$Day)
temp <- data.frame(temp)
#add those numbers to y data.frame; matched by Day
y$total <- temp$Freq[match(y$Day, temp$Var1)]
#create column colled percents in y that is frequency of cell / Day
y$percents <- (y$Freq/y$total)*100
#plot
g1 <- ggplot(data=subset(y, type %in% c('Fibro_1')), aes(x=Day, y=percents, group=type, color=type)) +
  geom_line(size=2)+
  geom_point(size=4) +
  xlab('Day post wounding') +
  ylab('% of all fibroblasts') +
  scale_color_manual(values=c('#F8766D')) +
  scale_y_continuous(limits = c(0, 60)) +
  theme_classic() +
  theme(legend.title = element_blank()) +
  theme(legend.position = c(0.7,0.85), legend.text = element_text(size=20)) +
  theme(axis.text.x = element_text(colour = "black", size = 20), axis.text.y = element_text(colour = 'black', size = 20)) +
  theme(axis.title.x = element_text(size=20), axis.title.y = element_text(size=20))
g2 <- ggplot(data=subset(y, type %in% c('Fibro_2')), aes(x=Day, y=percents, group=type, color=type)) +
  geom_line(size=2)+
  geom_point(size=4) +
  xlab('Day post wounding') +
  ylab('% of all fibroblasts') +
  scale_color_manual(values=c('#D89000')) +
  scale_y_continuous(limits = c(0, 60)) +
  theme_classic() +
  theme(legend.title = element_blank()) +
  theme(legend.position = c(0.7,0.85), legend.text = element_text(size=20)) +
  theme(axis.text.x = element_text(colour = "black", size = 20), axis.text.y = element_text(colour = 'black', size = 20)) +
  theme(axis.title.x = element_text(size=20), axis.title.y = element_text(size=20))
g3 <- ggplot(data=subset(y, type %in% c('Fibro_3')), aes(x=Day, y=percents, group=type, color=type)) +
  geom_line(size=2)+
  geom_point(size=4) +
  xlab('Day post wounding') +
  ylab('% of all fibroblasts') +
  scale_color_manual(values=c('#00BB4E')) +
  scale_y_continuous(limits = c(0, 60)) +
  theme_classic() +
  theme(legend.title = element_blank()) +
  theme(legend.position = c(0.7,0.85), legend.text = element_text(size=20)) +
  theme(axis.text.x = element_text(colour = "black", size = 20), axis.text.y = element_text(colour = 'black', size = 20)) +
  theme(axis.title.x = element_text(size=20), axis.title.y = element_text(size=20))
g4 <- ggplot(data=subset(y, type %in% c('Fibro_4')), aes(x=Day, y=percents, group=type, color=type)) +
  geom_line(size=2)+
  geom_point(size=4) +
  xlab('Day post wounding') +
  ylab('% of all fibroblasts') +
  scale_color_manual(values=c('#00BFC4')) +
  scale_y_continuous(limits = c(0, 60)) +
  theme_classic() +
  theme(legend.title = element_blank()) +
  theme(legend.position = c(0.7,0.85), legend.text = element_text(size=20)) +
  theme(axis.text.x = element_text(colour = "black", size = 20), axis.text.y = element_text(colour = 'black', size = 20)) +
  theme(axis.title.x = element_text(size=20), axis.title.y = element_text(size=20))
g5 <- ggplot(data=subset(y, type %in% c('Fibro_5')), aes(x=Day, y=percents, group=type, color=type)) +
  geom_line(size=2)+
  geom_point(size=4) +
  xlab('Day post wounding') +
  ylab('% of all fibroblasts') +
  scale_color_manual(values=c('magenta1')) +
  scale_y_continuous(limits = c(0, 60)) +
  theme_classic() +
  theme(legend.title = element_blank()) +
  theme(legend.position = c(0.5,0.85), legend.text = element_text(size=20)) +
  theme(axis.text.x = element_text(colour = "black", size = 20), axis.text.y = element_text(colour = 'black', size = 20)) +
  theme(axis.title.x = element_text(size=20), axis.title.y = element_text(size=20))
(g1 / g2 / g3 / g4 / g5)


#make Space/Time-tileplot Heatmaps -> Figure 3D (right)
#prep data.frame:
#generate data.frame 'y' with one column of all clusters and one column of all PunchID
y <- table(fibros@meta.data[,c('type', 'PunchID')])
y <- data.frame(y)
#make temporary data.frame with # of each cell per PunchID
temp <- table(fibros@meta.data$PunchID)
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
  p1<-ggplot(df, aes(x = Day, y=Space)) +
    geom_tile(aes(fill = percents_diff),color='black') +
    geom_text(aes(label = round(percents, 1)))+
    scale_fill_gradient2(low = 'blue',mid = 'white',high = 'red',midpoint = 0)+
    labs(title = paste(unique(blah$type)[i],'UW % :',as.character(round(uw_temp,1)),sep = ' '))
  plotlist[[i]]<-p1
}
grid.arrange(grobs = plotlist, ncol = 3) 

##### TilePlot cluster individually
df = subset(blah,subset = type == 'Fibro_1')
uw_temp = df$percents[df$Day=='UW']
df = subset(df,subset = Day != 'UW',invert=TRUE)
df$percents_diff = df$percents-uw_temp
tp1 <-ggplot(df, aes(x = Day, y=Space)) +
  geom_tile(aes(fill = percents_diff),color='black') +
  geom_text(aes(label = round(percents, 1)), size=7)+
  scale_fill_gradient2(low = 'blue',mid = 'white',high = 'red',midpoint = 0)+
  labs(title = paste('Fibro_1','UW% :',as.character(round(uw_temp,1)),sep = ' '))+
  coord_fixed()+
  theme_void()+  
  labs(fill='rel. to\nUW') +
  theme(axis.text = element_text(size = 20), plot.title = element_text(size=20), legend.title=element_text(size=20), legend.text = element_blank())
tp1
df = subset(blah,subset = type == 'Fibro_2')
uw_temp = df$percents[df$Day=='UW']
df = subset(df,subset = Day != 'UW',invert=TRUE)
df$percents_diff = df$percents-uw_temp
tp2 <-ggplot(df, aes(x = Day, y=Space)) +
  geom_tile(aes(fill = percents_diff),color='black') +
  geom_text(aes(label = round(percents, 1)), size=7)+
  scale_fill_gradient2(low = 'blue',mid = 'white',high = 'red',midpoint = 0)+
  labs(title = paste('fibros','UW% :',as.character(round(uw_temp,1)),sep = ' '))+
  coord_fixed()+
  theme_void()+
  labs(fill='rel. to\nUW') +
  theme(axis.text = element_text(size = 20), plot.title = element_text(size=20), legend.title=element_text(size=20), legend.text = element_blank())
tp2
df = subset(blah,subset = type == 'Fibro_3')
uw_temp = df$percents[df$Day=='UW']
df = subset(df,subset = Day != 'UW',invert=TRUE)
df$percents_diff = df$percents-uw_temp
tp3 <-ggplot(df, aes(x = Day, y=Space)) +
  geom_tile(aes(fill = percents_diff),color='black') +
  geom_text(aes(label = round(percents, 1)), size=7)+
  scale_fill_gradient2(low = 'blue',mid = 'white',high = 'red',midpoint = 0)+
  labs(title = paste('Fibro_3','UW% :',as.character(round(uw_temp,1)),sep = ' '))+
  coord_fixed()+
  theme_void()+
  labs(fill='rel. to\nUW') +
  theme(axis.text = element_text(size = 20), plot.title = element_text(size=20), legend.title=element_text(size=20), legend.text = element_blank())
tp3
df = subset(blah,subset = type == 'Fibro_4')
uw_temp = df$percents[df$Day=='UW']
df = subset(df,subset = Day != 'UW',invert=TRUE)
df$percents_diff = df$percents-uw_temp
tp4 <-ggplot(df, aes(x = Day, y=Space)) +
  geom_tile(aes(fill = percents_diff),color='black') +
  geom_text(aes(label = round(percents, 1)), size=7)+
  scale_fill_gradient2(low = 'blue',mid = 'white',high = 'red',midpoint = 0)+
  labs(title = paste('Fibro_4','UW% :',as.character(round(uw_temp,1)),sep = ' '))+
  coord_fixed()+
  theme_void()+
  labs(fill='rel. to\nUW') +
  theme(axis.text = element_text(size = 20), plot.title = element_text(size=20), legend.title=element_text(size=20), legend.text = element_blank())
tp4
df = subset(blah,subset = type == 'Fibro_5')
uw_temp = df$percents[df$Day=='UW']
df = subset(df,subset = Day != 'UW',invert=TRUE)
df$percents_diff = df$percents-uw_temp
tp5 <-ggplot(df, aes(x = Day, y=Space)) +
  geom_tile(aes(fill = percents_diff),color='black') +
  geom_text(aes(label = round(percents, 1)), size=7)+
  scale_fill_gradient2(low = 'blue',mid = 'white',high = 'red',midpoint = 0)+
  labs(title = paste('Fibro_5','UW% :',as.character(round(uw_temp,1)),sep = ' '))+
  coord_fixed()+
  theme_void()+
  labs(fill='rel. to\nUW') +
  theme(axis.text = element_text(size = 20), plot.title = element_text(size=20), legend.title=element_text(size=20), legend.text = element_blank())
tp5
pdf('fibro_tileplots.pdf', width = 4.4, height = 19.4)
(tp1 / tp2 / tp3 / tp4 / tp5)
dev.off()

#Jaccard similarity analysis -> Figures S3E and S3F
#######read in DEG lists of current study, Buechler et al. and Vu et al.
current.study.deg = readRDS("currentStudy.rds")
buechler.deg = readRDS("buechlerDEG.rds")
vu.deg = readRDS('vuDEG.rds')

####feed in lists of dataframes:
list1 = our.deg
list2 = buechler.deg
list3 = vu.deg

#generate matrix and plot Figure S3E: Current study vs Buechler et al.
jaccards = matrix(nrow = length(list1),ncol = length(list2))
###set cutoffs for p_val_adj & avg_logFC
for(n in 1:length(list1)){
  df1 = list1[[n]]
  df1$diff.pct = df1$pct.1 - df1$pct.2
  df1 = subset(df1,subset = p_val_adj<0.05)
  df1 = subset(df1,subset = avg_logFC>1) #the avg_logFC needs adjustment!
  #df1 = subset(df1,subset = `diff.pct`> 0.1)
  #df1 = df1[with(df1,order(-avg_logFC)),] #order by avg_logFC
  #df1 = df1[1:25,] #pick top 25 DEGs by previously ordered avg_logFC
  for(m in 1:length(list2)){
    #print(m)
    #print(df2$diff.pct)
    df2 = list2[[m]]
    df2$diff.pct = df2$pct.1 - df2$pct.2
    df2 = subset(df2,subset = p_val_adj<0.05)
    df2 = subset(df2,subset = avg_logFC>1) #the avg_logFC needs adjustment!
    #df2 = subset(df2,subset = `diff.pct`> 0.1)
    #df2 = df2[with(df2,order(-avg_log2FC)),] #order by avg_logFC
    #df2 = df2[1:25,] #pick top 25 DEGs by previously ordered avg_logFC
    #beepr::beep(2)
    
    #calculate jaccards
    j.temp = length(intersect(df1$gene,df2$Gene)) / length(union(df1$gene,df2$Gene)) #careful about upper/lower case 'Gene'!
    jaccards[n,m] = j.temp
    
  }
}
row.names(jaccards) = names(list1)
colnames(jaccards) = names(list2)
jaccards
#create heatmap, first melt matrix
library(reshape2)
melt_jaccards <- melt(jaccards)
head(melt_jaccards)
melt_jaccards$Var1 = as.character(melt_jaccards$Var1)
melt_jaccards$Var2 = as.character(melt_jaccards$Var2)

#create heatmap
pdf('jaccard_buechler.pdf', width = 5, height = 4.5)
ggplot(data = melt_jaccards, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  scale_fill_gradient2(low = 'blue', high='red', mid = 'white', midpoint = 0) +
  theme_classic() +
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text = element_text(colour = 'black'))
dev.off()

#generate matrix and plot Figure S3F: Current study vs Vu et al.
jaccards = matrix(nrow = length(list1),ncol = length(list3))
###set cutoffs for p_val_adj & avg_logFC
for(n in 1:length(list1)){
  df1 = list1[[n]]
  df1$diff.pct = df1$pct.1 - df1$pct.2
  df1 = subset(df1,subset = p_val_adj<0.05)
  df1 = subset(df1,subset = avg_logFC>1) #the avg_logFC needs adjustment!
  #df1 = subset(df1,subset = `diff.pct`> 0.1)
  #df1 = df1[with(df1,order(-avg_logFC)),] #order by avg_logFC
  #df1 = df1[1:25,] #pick top 25 DEGs by previously ordered avg_logFC
  for(m in 1:length(list3)){
    #print(m)
    #print(df2$diff.pct)
    df2 = list3[[m]]
    df2$diff.pct = df2$pct.1 - df2$pct.2
    df2 = subset(df2,subset = p_val_adj<0.05)
    df2 = subset(df2,subset = avg_logFC>1) #the avg_logFC needs adjustment!
    #df2 = subset(df2,subset = `diff.pct`> 0.1)
    #df2 = df2[with(df2,order(-avg_log2FC)),] #order by avg_logFC
    #df2 = df2[1:25,] #pick top 25 DEGs by previously ordered avg_logFC
    #beepr::beep(2)
    
    #calculate jaccards
    j.temp = length(intersect(df1$gene,df2$gene)) / length(union(df1$gene,df2$gene)) #careful about upper/lower case 'Gene'!
    jaccards[n,m] = j.temp
    
  }
}
row.names(jaccards) = names(list1)
colnames(jaccards) = names(list3)
jaccards
#create heatmap, first melt matrix
library(reshape2)
melt_jaccards <- melt(jaccards)
head(melt_jaccards)
melt_jaccards$Var1 = as.character(melt_jaccards$Var1)
melt_jaccards$Var2 = as.character(melt_jaccards$Var2)

#create heatmap
pdf('jaccard_vu.pdf', width = 5, height = 4.5)
ggplot(data = melt_jaccards, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  scale_fill_gradient2(low = 'blue', high='red', mid = 'white', midpoint = 0) +
  theme_classic() +
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text = element_text(colour = 'black'))
dev.off()