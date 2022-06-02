
####################################################################
#On to correlation matrix calculation: ##########
################################################
#load in all the stuff
mm = readRDS("E:/Krummel Lab/SC_Experiments/20201130_WH011/WH_REDO/nmf/combinations/nmf_all_cells/mm_24_ns.rds")
neuts = readRDS("D:/Krummel Lab/SC_Experiments/20201130_WH011/WH_REDO/nmf/combinations/nmf_all_cells/neut_9_s.rds")
t = readRDS("D:/Krummel Lab/SC_Experiments/20201130_WH011/WH_REDO/nmf/combinations/nmf_all_cells/t_10_s.rds")
tnk = readRDS("D:/Krummel Lab/SC_Experiments/20201130_WH011/WH_REDO/nmf/combinations/nmf_all_cells/tnk_3_s.rds")
dc = readRDS("D:/Krummel Lab/SC_Experiments/20201130_WH011/WH_REDO/nmf/combinations/nmf_all_cells/dc_6_s.rds")
b = readRDS("D:/Krummel Lab/SC_Experiments/20201130_WH011/WH_REDO/nmf/combinations/nmf_all_cells/b_2_s.rds")
mast = readRDS("D:/Krummel Lab/SC_Experiments/20201130_WH011/WH_REDO/nmf/combinations/nmf_all_cells/mast_2_s.rds")
fibro = readRDS("D:/Krummel Lab/SC_Experiments/20201130_WH011/WH_REDO/nmf/combinations/nmf_all_cells/fibro_17_ns.rds")
endo = readRDS("D:/Krummel Lab/SC_Experiments/20201130_WH011/WH_REDO/nmf/combinations/nmf_all_cells/endo_11_ns.rds")
kerat = readRDS("D:/Krummel Lab/SC_Experiments/20201130_WH011/WH_REDO/nmf/combinations/nmf_all_cells/kerat_16_s.rds")
vsm = readRDS("D:/Krummel Lab/SC_Experiments/20201130_WH011/WH_REDO/nmf/combinations/nmf_all_cells/vsm_3_s.rds")
dsp = readRDS("D:/Krummel Lab/SC_Experiments/20201130_WH011/WH_REDO/nmf/combinations/nmf_all_cells/dsp_6_s.rds")
melano = readRDS("D:/Krummel Lab/SC_Experiments/20201130_WH011/WH_REDO/nmf/combinations/nmf_all_cells/melano_5_s.rds")
test = rbind(mm,neuts,t,tnk,dc,b,mast,fibro,endo,kerat,vsm,dsp,melano)
########3
library(scales)
library(reshape2)
blah = readRDS("E:/Krummel Lab/SC_Experiments/20201130_WH011/WH_REDO/nmf/combinations/mm_nmf_scale_table.rds")
blah2 = readRDS("E:/Krummel Lab/SC_Experiments/20201130_WH011/WH_REDO/nmf/combinations/fibro_nmf_scale_table.rds")
blah = readRDS("E:/Krummel Lab/SC_Experiments/20201130_WH011/WH_REDO/nmf/combinations/mm_nmf_table.rds")
blah2 = readRDS("E:/Krummel Lab/SC_Experiments/20201130_WH011/WH_REDO/nmf/combinations/cd45pos_freqs_table.rds")
blah2 = readRDS("E:/Krummel Lab/SC_Experiments/20201130_WH011/WH_REDO/nmf/combinations/fibro_nmf_table.rds")
blah = mm
rownames(blah) = paste('MM ',rownames(blah))
blah2 = blah.pos
blah = blah[,-c(13)]
blah2 = blah2[,-c(13)]

big.corr_cross = cor(t(blah),t(blah2),method = 'pearson')
d = dist(big.corr_cross)
hc1 = hclust(d)
variables1_reorder = hc1$labels[hc1$order]

d = dist(t(big.corr_cross))
hc2 = hclust(d)
variables2_reorder = hc2$labels[hc2$order]

corr_df <-as.data.frame(big.corr_cross)
colnames(corr_df)<-hc2$labels
corr_df = corr_df[,variables2_reorder]
corr_df$id = factor(hc1$labels,levels = variables1_reorder)

melted = melt(corr_df,id.vars = "id")
melted = melted[order(melted$id),]
library(ggplot2)
library(scales)
ggplot(data = melted,aes(x=id, y=variable, fill=value))+geom_tile()+
  scale_fill_gradient2(low = muted("blue"),high = muted("red"),mid = "white",midpoint=0,na.value = "grey20",limits = c(-1,1))+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.4),axis.title.x = element_blank(),axis.title.y = element_blank())+coord_fixed()
##########for significance testing:##############
###############################################
melted$label = rep(" ",nrow(melted))
melted$pval = integer(nrow(melted))+1
melted$pval_adj = integer(nrow(melted))+1
for (i in 1:nrow(melted)){
  derp<-cor.test(blah2[match(as.character(melted$variable[i]),row.names(blah2)),] %>% as.numeric(),
                 blah[match(as.character(melted$id[i]),row.names(blah)),] %>% as.numeric(),method = "pearson")
  melted$pval[i]<-derp$p.value
}

melted$label[melted$pval < (0.05)] = '+'
melted$label[melted$pval < (0.005)] = '++'
ggplot(data = melted,aes(x=id, y=variable, fill=value))+geom_tile()+
  scale_fill_gradient2(low = muted("blue"),high = muted("red"),mid = "white",midpoint=0,na.value = "grey20",limits = c(-1,1))+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.4),axis.title.x = element_blank(),axis.title.y = element_blank())+geom_text(aes(label=label), color="black", size=3,nudge_y = 0)+
  coord_fixed()

####################################################
melted$label_adj = rep(" ",nrow(melted))
temp = p.adjust(melted$pval, method = "BH")
melted$pval_adj=temp
melted$label_adj[melted$pval_adj < (0.05)] = '+'
melted$label_adj[melted$pval_adj < (0.005)] = '++'
ggplot(data = melted,aes(x=id, y=variable, fill=value))+geom_tile()+
  scale_fill_gradient2(low = muted("blue"),high = muted("red"),mid = "white",midpoint=0,na.value = "grey20",limits = c(-1,1))+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.4,size=12),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.y = element_text(size=12))+
  geom_text(aes(label=label_adj), color="black", size=4,nudge_y = 0)+
  coord_fixed()
###################################################
irds = (table(cd45.neg.sub$PunchID))
cd45.neg.sub$Punch_ID <- factor(cd45.neg.sub$PunchID, levels = names(irds))
Idents(cd45.neg.sub)<-cd45.neg.sub$Punch_ID
cluster.averages <- AverageExpression(cd45.neg.sub,assays = 'NMF',slot = 'counts' )
blah2 = cluster.averages$NMF
#######
big.corr_cross = cor(t(blah),t(blah2),method = 'pearson')
d = dist(big.corr_cross)
hc1 = hclust(d)
variables1_reorder = hc1$labels[hc1$order]

d = dist(t(big.corr_cross))
hc2 = hclust(d)
variables2_reorder = hc2$labels[hc2$order]

corr_df <-as.data.frame(big.corr_cross)
colnames(corr_df)<-hc2$labels
corr_df = corr_df[,variables2_reorder]
corr_df$id = factor(hc1$labels,levels = variables1_reorder)

melted = melt(corr_df,id.vars = "id")
melted = melted[order(melted$id),]

ggplot(data = melted,aes(x=id, y=variable, fill=value))+geom_tile()+
  scale_fill_gradient2(low = muted("blue"),high = muted("red"),mid = "white",midpoint=0,na.value = "grey20",limits = c(-1,1))+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.4),axis.title.x = element_blank(),axis.title.y = element_blank())+coord_fixed()

melted$label = rep(" ",nrow(melted))
melted$pval = integer(nrow(melted))+1
melted$pval_adj = integer(nrow(melted))+1
for (i in 1:nrow(melted)){
  derp<-cor.test(blah2[match(as.character(melted$variable[i]),row.names(blah2)),] %>% as.numeric(),
                 blah[match(as.character(melted$id[i]),row.names(blah)),] %>% as.numeric(),method = "pearson")
  melted$pval[i]<-derp$p.value
  
  
}
temp = p.adjust(melted$pval, method = "BH")
melted$pval_adj=temp
melted$label[melted$pval < (0.05)] = '+'
melted$label[melted$pval < (0.005)] = '++'
ggplot(data = melted,aes(x=id, y=variable, fill=value))+geom_tile()+
  scale_fill_gradient2(low = muted("blue"),high = muted("red"),mid = "white",midpoint=0,na.value = "grey20",limits = c(-1,1))+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.4),axis.title.x = element_blank(),axis.title.y = element_blank())+geom_text(aes(label=label), color="black", size=3,nudge_y = 0)+
  coord_fixed()
corr = matrix(0,ncol(w),ncol(tumor_w))
row.names(corr) = paste('factor-',c(1:ncol(wh_w)),sep='')
colnames(corr) =  paste('factor-',c(1:ncol(tumor_w)),sep='')
for(i in 1:ncol(wh_w)){
  temp_wh = wh_w[,i]
  temp_wh = temp_wh[match(mutual_genes,names(temp_wh))]
  for(j in 1:ncol(tumor_w)){
    temp_tumor = tumor_w[,j]
    temp_tumor = temp_tumor[match(mutual_genes,names(temp_tumor))]
    corr[i,j] = cor(temp_wh,temp_tumor,method = 'pearson')
  }
}
####for big correlation for all factors####
###fo####################################33333333
##combine into bigass correlation matrix symmetrical:

big.blah = test
big.corr = cor(t(big.blah),method = "pearson")
#big.corr = cor(t(big.blah),method = "spearman")
d = dist(big.corr)
hc1 = hclust(d)
variables1_reorder = hc1$labels[hc1$order]
corr_df <-as.data.frame(big.corr)
colnames(corr_df)<-hc1$labels
corr_df = corr_df[,variables1_reorder]
corr_df$id = factor(hc1$labels,levels = variables1_reorder)
melted = melt(corr_df,id.vars = "id")
melted = melted[order(melted$id),]

ggplot(data = melted,aes(x=id, y=variable, fill=value))+geom_tile()+
  scale_fill_gradient2(low = muted("blue"),high = muted("red"),mid = "white",midpoint=0,na.value = "grey20",limits = c(-1,1))+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.4),axis.title.x = element_blank(),axis.title.y = element_blank())+
  coord_fixed()

####try with corrplot, Figure S6H
library(corrplot)
#calculate rhos
big.blah = test

c = cor(t(big.blah),method = "pearson")
#calculate pval matrix
p=c
for(k in 1:nrow(c)){
  for(j in 1:ncol(c)){
    p[k,j] = cor.test(as.matrix(test)[k,],as.matrix(test)[j,],method = "pearson")$p.value
  }
}
#trim labels
row.names(c)=  str_sub(row.names(c),3,-1)
colnames(c)=  str_sub(colnames(c),3,-1)
row.names(p)=  str_sub(row.names(p),3,-1)
colnames(p)=  str_sub(colnames(p),3,-1)
cp = corrplot(c, p.mat = p, sig.level = 0.05, insig = "blank",col = rev(brewer.pal(10 ,"RdBu")),tl.col = 'black',order = 'hclust',addrect = 25)

###normalize to max value across rows
temp.blah = big.blah
for (i in 1:nrow(big.blah)){
  temp.blah[i,] = temp.blah[i,]/max(temp.blah[i,])
}
temp.blah$celltype = gsub( " .*$", "", row.names(temp.blah) )
temp.blah$factor = factor(row.names(temp.blah),levels = row.names(cp$corr))
ggplot(temp.blah,aes(x=factor,y=D01_2mm))+geom_bar(stat = 'identity',aes(fill=celltype))
temp.blah = temp.blah.master
#######
#make interpolations
d13 = matrix(nrow = nrow(temp.blah),ncol = 8)
row.names(d13) = row.names(temp.blah)
d13[,c(1:4)] = ((as.matrix(temp.blah[,c(5:8)])-as.matrix(temp.blah[,c(1:4)]))/3+as.matrix(temp.blah[,c(1:4)])) %>% as.numeric()
d13[,c(5:8)] = ((as.matrix(temp.blah[,c(5:8)])-as.matrix(temp.blah[,c(1:4)]))/3*2+as.matrix(temp.blah[,c(1:4)])) %>% as.numeric()

d37 = matrix(nrow = nrow(temp.blah),ncol = 8)
row.names(d37) = row.names(temp.blah)
d37[,c(1:4)] = ((as.matrix(temp.blah[,c(9:12)])-as.matrix(temp.blah[,c(5:8)]))/3+as.matrix(temp.blah[,c(5:8)])) %>% as.numeric()
d37[,c(5:8)] = ((as.matrix(temp.blah[,c(9:12)])-as.matrix(temp.blah[,c(5:8)]))/3*2+as.matrix(temp.blah[,c(5:8)])) %>% as.numeric()

d714 = matrix(nrow = nrow(temp.blah),ncol = 8)
row.names(d13) = row.names(temp.blah)
d714[,c(1:4)] = ((as.matrix(temp.blah[,c(13:16)])-as.matrix(temp.blah[,c(9:12)]))/3+as.matrix(temp.blah[,c(9:12)])) %>% as.numeric()
d714[,c(5:8)] = ((as.matrix(temp.blah[,c(13:16)])-as.matrix(temp.blah[,c(9:12)]))/3*2+as.matrix(temp.blah[,c(9:12)])) %>% as.numeric()

bold.lines = c(0,10,16,20,24,25,30,39,40,45,47,52,61,65,67,71,78,81,88,90,98,105,111,112,113,114)
library(scales)
cols.use = (hue_pal()(13))
circos.clear()
circos.par("track.height" = 0.125,"gap.degree"=0,'start.degree'=90,'cell.padding'=c(0.0, 1.00, 0.0, 1.00))
circos.initialize(letters[1:1], xlim = c(0, 114))
#outer track, blank 
circos.track(ylim = c(0, 1), bg.border = 'white',panel.fun = function(x, y) {
  value = rep(0,114)
  circos.barplot(value, border='white',pos = as.numeric(temp.blah$factor)-0.5,bar_width=1,col = cols.use[as.factor(temp.blah$celltype)])
})
set_track_gap(0) 
circos.axis("bottom",major.at =bold.lines,minor.ticks = 0,major.tick.length = 1,labels.cex = 0,lwd=3)

#UW track
circos.track(ylim = c(0, 1), track.height = 0.08,panel.fun = function(x, y) {
  value = temp.blah.master[,17]
  circos.barplot(value, pos = as.numeric(temp.blah$factor)-0.5,bar_width=1,col = cols.use[as.factor(temp.blah$celltype)])
})
#add labels
circos.text(x=1:114-0.5,y=rep(1.2,114),cex=0.8 ,facing = "clockwise",niceFacing = TRUE,labels = row.names(cp$corr)[1:114],adj = c(0,0.5),)
#circos.text(x=58:114-1.2,y=rep(1,57), labels = row.names(cp$corr)[58:114],adj = c(0,-1),facing = "reverse.clockwise")
circos.axis("bottom",major.at =0:114,minor.ticks = 0,major.tick.length = 1,labels.cex = 0,lwd=0.5)
circos.axis("bottom",major.at =bold.lines,minor.ticks = 0,major.tick.length = 1,labels.cex = 0,lwd=3)
set_track_gap(0.03) 
#outermost ring
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  value = temp.blah[,16]
  circos.barplot(value, pos = as.numeric(temp.blah$factor)-0.5,bar_width=1,col = cols.use[as.factor(temp.blah$celltype)])
})
circos.axis("bottom",major.at =0:114,minor.ticks = 0,major.tick.length = 1,labels.cex = 0,lwd=0.5)
circos.axis("bottom",major.at =bold.lines,minor.ticks = 0,major.tick.length = 1,labels.cex = 0,lwd=3)

set_track_gap(0.02) 
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  value = temp.blah[,15]
  circos.barplot(value, pos = as.numeric(temp.blah$factor)-0.5,bar_width=1,col = cols.use[as.factor(temp.blah$celltype)])
})
circos.axis("bottom",major.at =1:114,minor.ticks = 0,major.tick.length = 1,labels.cex = 0,lwd=0.5)
circos.axis("bottom",major.at =bold.lines,minor.ticks = 0,major.tick.length = 1,labels.cex = 0,lwd=3)

circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  value = temp.blah[,14]
  circos.barplot(value, pos = as.numeric(temp.blah$factor)-0.5,bar_width=1,col = cols.use[as.factor(temp.blah$celltype)])
})
circos.axis("bottom",major.at =1:114,minor.ticks = 0,major.tick.length = 1,labels.cex = 0,lwd=0.5)
circos.axis("bottom",major.at =bold.lines,minor.ticks = 0,major.tick.length = 1,labels.cex = 0,lwd=3)

circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  value = temp.blah[,13]
  circos.barplot(value, pos = as.numeric(temp.blah$factor)-0.5,bar_width=1,col = cols.use[as.factor(temp.blah$celltype)])
})
circos.axis("bottom",major.at =1:114,minor.ticks = 0,major.tick.length = 1,labels.cex = 0,lwd=0.5)
circos.axis("bottom",major.at =bold.lines,minor.ticks = 0,major.tick.length = 1,labels.cex = 0,lwd=3)

###############################################3
#for interpolated ones
interpol = d714
bold.lines = c(0,10,16,20,24,25,30,39,40,45,47,52,61,65,67,71,78,81,88,90,98,105,111,112,113,114)
library(scales)
cols.use = (hue_pal()(13))
circos.clear()
circos.par("track.height" = 0.125,"gap.degree"=0,'start.degree'=90,'cell.padding'=c(0.0, 1.00, 0.0, 1.00))
circos.initialize(letters[1:1], xlim = c(0, 114))
#outer track, blank 
circos.track(ylim = c(0, 1), bg.border = 'white',panel.fun = function(x, y) {
  value = rep(0,114)
  circos.barplot(value, border='white',pos = as.numeric(temp.blah$factor)-0.5,bar_width=1,col = cols.use[as.factor(temp.blah$celltype)])
})
set_track_gap(0) 
circos.axis("bottom",major.at =bold.lines,minor.ticks = 0,major.tick.length = 1,labels.cex = 0,lwd=3)

#UW track
circos.track(ylim = c(0, 1), track.height = 0.08,panel.fun = function(x, y) {
  value = temp.blah[,17]
  circos.barplot(value, pos = as.numeric(temp.blah$factor)-0.5,bar_width=1,col = cols.use[as.factor(temp.blah$celltype)])
})
#add labels
circos.text(x=1:114-0.5,y=rep(1.2,114),cex=0.8 ,facing = "clockwise",niceFacing = TRUE,labels = row.names(cp$corr)[1:114],adj = c(0,0.5),)
#circos.text(x=58:114-1.2,y=rep(1,57), labels = row.names(cp$corr)[58:114],adj = c(0,-1),facing = "reverse.clockwise")
circos.axis("bottom",major.at =0:114,minor.ticks = 0,major.tick.length = 1,labels.cex = 0,lwd=0.5)
circos.axis("bottom",major.at =bold.lines,minor.ticks = 0,major.tick.length = 1,labels.cex = 0,lwd=3)
set_track_gap(0.03) 
#outermost ring
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  value = interpol[,8]
  circos.barplot(value, pos = as.numeric(temp.blah$factor)-0.5,bar_width=1,col = cols.use[as.factor(temp.blah$celltype)])
})
circos.axis("bottom",major.at =0:114,minor.ticks = 0,major.tick.length = 1,labels.cex = 0,lwd=0.5)
circos.axis("bottom",major.at =bold.lines,minor.ticks = 0,major.tick.length = 1,labels.cex = 0,lwd=3)

set_track_gap(0.02) 
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  value = interpol[,7]
  circos.barplot(value, pos = as.numeric(temp.blah$factor)-0.5,bar_width=1,col = cols.use[as.factor(temp.blah$celltype)])
})
circos.axis("bottom",major.at =1:114,minor.ticks = 0,major.tick.length = 1,labels.cex = 0,lwd=0.5)
circos.axis("bottom",major.at =bold.lines,minor.ticks = 0,major.tick.length = 1,labels.cex = 0,lwd=3)

circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  value = interpol[,6]
  circos.barplot(value, pos = as.numeric(temp.blah$factor)-0.5,bar_width=1,col = cols.use[as.factor(temp.blah$celltype)])
})
circos.axis("bottom",major.at =1:114,minor.ticks = 0,major.tick.length = 1,labels.cex = 0,lwd=0.5)
circos.axis("bottom",major.at =bold.lines,minor.ticks = 0,major.tick.length = 1,labels.cex = 0,lwd=3)

circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  value = interpol[,5]
  circos.barplot(value, pos = as.numeric(temp.blah$factor)-0.5,bar_width=1,col = cols.use[as.factor(temp.blah$celltype)])
})
circos.axis("bottom",major.at =1:114,minor.ticks = 0,major.tick.length = 1,labels.cex = 0,lwd=0.5)
circos.axis("bottom",major.at =bold.lines,minor.ticks = 0,major.tick.length = 1,labels.cex = 0,lwd=3)


###########PLOT ALL 114 TILEPLOTS for FIgure S6I ####
######################TILEPLOT FOR ALL FACTORS$##################
#obj = cd45.pos.sub
library(gridExtra)
library(ggplot2)
library(magrittr)
library(reshape2)

plotlist = list()
factor.order = cp$corr %>% row.names()
for (i in 1:length(factor.order)){
  feature.interest = row.names(big.blah)[match(factor.order[i],row.names(big.blah))]
  print(match(factor.order[i],row.names(big.blah)))
  temp = big.blah[match(factor.order[i],row.names(big.blah)),]
  feature.averages = melt(temp)
  feature.averages$Day = sub("\\_.*", "", feature.averages$variable)
  feature.averages$Space = sub(".*\\_","",feature.averages$variable)
  
  uw_temp = feature.averages$value[match('UW',feature.averages$variable)]
  df = subset(feature.averages,subset = Day != 'UW',invert=TRUE)
  df$diff = df$value-uw_temp
  rm(p1)
  p1<-ggplot(df, aes(x = Day, y=Space)) +
    geom_tile(aes(fill = diff),color='black') +
    geom_text(aes(label = round(value, 1)),size=4)+
    scale_fill_gradient2(low = 'blue',mid = 'white',high = 'red',midpoint = 0)+
    labs(title = paste(feature.interest,'UW:',as.character(round(uw_temp,1)),sep = ' '))+coord_fixed()+
    theme_void()+theme(axis.text = element_text(size = 10),legend.title = element_blank(),legend.position = 'None',
                       plot.margin = margin(t=20,r=20,b=20,l=20))
  
  plotlist[[i]]<-p1
  
}
grid.arrange(grobs = plotlist, ncol = 12) 
################################################################