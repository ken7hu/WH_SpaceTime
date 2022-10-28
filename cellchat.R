#script for Cellchat analyssis of interactions between the mono.mac subset and fibro subsets
#start with cd45.full from the completely merged dataset-prolifs-mitohigh macs 
cells.desired = c('Fibro_1','Fibro_2','Fibro_3','Fibro_4',
                  'Fibro_5','Mono','Mono_Mac_1','Mono_Mac_2','Mono_Mac_4','Mono_Mac_5','Mono_Mac_6',
                  'Mono_Mac_Mgl2','Mono_Mac_MHCII','Mono_Vcan')
######redo, merge the fibro and mm objects
cd45.full <- merge(cd45.pos.sub, y = cd45.neg.sub, add.cell.ids = c("MM", "FIBRO"), project = "WH")

all.subset = subset(cd45.full,idents = cells.desired)
#here we go...
test_list = SplitObject(cd45.full,split.by = 'Day')
#1 is d01
#2 is d07
#3 is d00
#4 is d03
#5 is d14
DB = CellChatDB.mouse

cellchat_01 = createCellChat(object = test_list[[1]],group.by = 'ident',assay = 'RNA')
cellchat_01@DB <- DB
cellchat_01 <- subsetData(cellchat_01)
future::plan("multiprocess",workers=6)
cellchat_01<- identifyOverExpressedGenes(cellchat_01)
cellchat_01 <- identifyOverExpressedInteractions(cellchat_01)
beepr::beep(2)
cellchat_01<-computeCommunProb(cellchat_01,type = "triMean",raw.use = TRUE,population.size = TRUE)
cellchat_01 <- computeCommunProbPathway(cellchat_01)
cellchat_01 <- aggregateNet(cellchat_01)
cellchat_01 <- netAnalysis_computeCentrality(cellchat_01,slot.name = "netP")
###
cellchat_07 = createCellChat(object = test_list[[2]],group.by = 'ident',assay = 'RNA')
cellchat_07@DB <- DB
cellchat_07 <- subsetData(cellchat_07)
future::plan("multiprocess",workers=6)
cellchat_07 <- identifyOverExpressedGenes(cellchat_07)
cellchat_07 <- identifyOverExpressedInteractions(cellchat_07)
beepr::beep(2)
cellchat_07<-computeCommunProb(cellchat_07,type = "triMean",raw.use = TRUE,population.size = TRUE)
cellchat_07 <- computeCommunProbPathway(cellchat_07)
cellchat_07 <- aggregateNet(cellchat_07)
cellchat_07 <- netAnalysis_computeCentrality(cellchat_07,slot.name = "netP")
#####adding in cellchat 03
cellchat_03 = createCellChat(object = test_list[[4]],group.by = 'ident',assay = 'RNA')
cellchat_03@DB <- DB
cellchat_03 <- subsetData(cellchat_03)
future::plan("multiprocess",workers=6)
cellchat_03<- identifyOverExpressedGenes(cellchat_03)
cellchat_03 <- identifyOverExpressedInteractions(cellchat_03)
cellchat_03<-computeCommunProb(cellchat_03,type = "triMean",raw.use = TRUE,population.size = TRUE)
cellchat_03 <- computeCommunProbPathway(cellchat_03)
cellchat_03 <- aggregateNet(cellchat_03)
cellchat_03 <- netAnalysis_computeCentrality(cellchat_03,slot.name = "netP")
########d14
cellchat_14 = createCellChat(object = test_list[[5]],group.by = 'ident',assay = 'RNA')
cellchat_14@DB <- DB
cellchat_14 <- subsetData(cellchat_14)
future::plan("multiprocess",workers=6)
cellchat_14 <- identifyOverExpressedGenes(cellchat_14)
cellchat_14 <- identifyOverExpressedInteractions(cellchat_14)
cellchat_14<-computeCommunProb(cellchat_14,type = "triMean",raw.use = TRUE,population.size = TRUE)
cellchat_14 <- computeCommunProbPathway(cellchat_14)
cellchat_14 <- aggregateNet(cellchat_14)
cellchat_14 <- netAnalysis_computeCentrality(cellchat_14,slot.name = "netP")
######for d00
cellchat_00 = createCellChat(object = test_list[[3]],group.by = 'ident',assay = 'RNA')
cellchat_00@DB <- DB
cellchat_00 <- subsetData(cellchat_00)
future::plan("multiprocess",workers=6)
cellchat_00<- identifyOverExpressedGenes(cellchat_00)
cellchat_00 <- identifyOverExpressedInteractions(cellchat_00)
beepr::beep(2)
cellchat_00<-computeCommunProb(cellchat_00,type = "triMean",raw.use = TRUE,population.size = TRUE)
cellchat_00 <- computeCommunProbPathway(cellchat_00)
cellchat_00 <- aggregateNet(cellchat_00)
cellchat_00 <- netAnalysis_computeCentrality(cellchat_00,slot.name = "netP")

#######################
object.list <- list(D01 = cellchat_01,D03 = cellchat_03,D07 = cellchat_07,D14 = cellchat_14,D00 = cellchat_00)
saveRDS(object.list,file = "E:/Krummel Lab/SC_Experiments/20201130_WH011/WH_REDO/cellchat/mac_and_fibro_only/20210909_mac_fibro_cellchatlist_reannotated.rds")
#saved this under 20210531cellchatobjectlist
cellchat_big <- mergeCellChat(object.list, add.names = names(object.list))
netVisual_heatmap(cellchat_big, comparison = c(1,4) )
netVisual_heatmap(cellchat_big,measure = "weight",comparison = c(1,5))
##################
#stacked barchart## (Supplementary Figure 5D)
##################
df2<-rankNet(cellchat_big, mode = "comparison",stacked = T,do.stat=FALSE,comparison = c(1,2,3,4,5),return.data = TRUE)
temp2 = df2$signaling.contribution
names.list = unique(temp2$name)
weight.avg = c()
for (i in 1:length(names.list)){
  weight.avg[i]=sum(temp2[temp2$name == names.list[i],]$contribution * c(1,2,3,4,5)) / (sum(temp2[temp2$name == names.list[i],]$contribution))
}
temp2$name = factor(temp2$name,levels = temp2$name[order(weight.avg)])
temp2$group = factor(temp2$group,levels = c('D00','D14','D07','D03','D01'))

gg <- ggplot(temp2, aes(x=name, y=contribution, fill = group)) + geom_bar(stat="identity",width = 0.75, position ="fill") +
  xlab("") + ylab("Relative information flow") + coord_flip()#+ theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
#  scale_y_discrete(breaks=c("0","0.5","1")) +
gg <- gg + geom_hline(yintercept = 0.5, linetype="dashed", color = "grey50", size=0.5)+scale_y_continuous(expand = c(0, 0))+theme_classic()+
  theme(legend.text = element_text(size=16),legend.title = element_blank())
gg

##############3# Hierarchy plot Figure S5E
pathways.show <- c("PERIOSTIN") 
weight.max = 0.00065531/4
vertex.receiver = seq(1,5) # Left portion of hierarchy plot the shows signaling to dermal cells and right portion shows signaling to epidermal cells
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:c(2)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, vertex.receiver = vertex.receiver, edge.weight.max = weight.max, edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
for (i in 1:c(3)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, vertex.receiver = vertex.receiver, edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
