############################################################################
######feed in nmf parameter sweep for number of factors
nmf_outs = list()
indices = c(18:26)
cd45.pos.sub = readRDS("E:/Krummel Lab/SC_Experiments/20201130_WH011/WH_REDO/nmf/cd45neg_nmf/")
cd45.pos.sub = readRDS("E:/Krummel Lab/SC_Experiments/20201130_WH011/WH_REDO/nmf/cd45pos/second_try/neuts_nmf/20210930_neuts_subset.rds")
for (i in c(1:length(indices))){
  print(i)
  nmf_outs[[i]] = readRDS(paste("~/Krummel lab/20201130_WH011/WH_REDO/nmf/humu/nmf_seurat_scaling_wh_ref_genes/nmf_out_",
                                as.character(indices[i]),'.rds',sep=''))
  
}
for (i in c(1:length(indices))){
  print(i)
  nmf_outs[[i]] = readRDS(paste("~/Krummel lab/20201130_WH011/WH_REDO/nmf/cd45pos/second_try/mm_outs_seurat_scaling/nmf_ns_",
                                as.character(indices[i]),'.rds',sep=''))
}
names(nmf_outs) = paste(as.character(indices),'_factors',sep='')
#iterate through nnf_outs to calculate cophcor and dispersion
cophenetic = c()
dispersion = c()
for(j in 1:length(nmf_outs)){
  print(j)
  cophenetic[j] = cophcor(nmf_outs[[j]])
  dispersion[j] = dispersion(nmf_outs[[j]])
}
df = data.frame(indices,cophenetic,dispersion)
####plotting metrics cophenetic and dispersion score
df = df[3:11,]
colnames(df)[1]<-'factors'
ggplot(data = df, aes(x=factors,y=cophenetic))+geom_line(size=1)+geom_point(size=3)+theme_bw()+
  theme(axis.text = element_text(size=12),axis.title = element_text(size=12))
########
par(mfrow = c(2,1), xpd=TRUE)
plot(x=indices,y=cophenetic)
plot(x=indices,y=dispersion)
#########################################################
##### for analyzing selected NMF output
nmf.out.24 = readRDS("~/Krummel lab/20201130_WH011/WH_REDO/nmf/cd45pos/second_try/mm_outs_seurat_scaling/nmf_ns_24.rds")
w<-basis(nmf.out.24)
h<-coef(nmf.out.24)
row.names(h)<-paste('factor-',as.character(c(1:nrow(h))),sep ='')
cd45.pos.sub[["NMF"]] <- CreateAssayObject(counts = h)
FeaturePlot(cd45.pos.sub,features = c('factor-1','factor-2','factor-3','factor-4','factor-5','factor-6','factor-7',"factor-8","factor-9","factor-10",
                                      "factor-11","factor-12","factor-13","factor-14","factor-15","factor-16"),ncol = 4,slot = 'counts',reduction = 'umap')
FeaturePlot(cd45.pos.sub,features = c('factor-17','factor-18','factor-19','factor-20','factor-21','factor-22',"factor-23","factor-24","factor-25",
                                      "factor-26","factor-27","factor-28","factor-29","factor-30"),ncol = 4,slot = 'counts',reduction = 'umap')

#####PLOT ALL THE TILE PLOTS FOR A subtype: ##########
#####################################################
library(gridExtra)
plotlist = list()
poo = (table(obj$PunchID))
obj$Punch_ID <- factor(obj$PunchID, levels = names(poo))
Idents(obj)<-obj$Punch_ID
for (i in 1:nrow(obj@assays$NMF)){
  print(i)
  feature.interest = row.names(obj@assays$NMF)[i]
  feature.averages = AverageExpression(obj,assays = 'NMF',slot = 'counts',features = feature.interest,verbose = 'FALSE')
  #feature.averages = aggregate(obj[[feature]], list(obj$Punch_ID), mean)
  feature.averages = as.data.frame(t(feature.averages[[1]]))
  feature.averages$Day = sub("\\_.*", "", row.names(feature.averages))
  feature.averages$Space = sub(".*\\_","",row.names(feature.averages))
  uw_temp = feature.averages[[feature.interest]][feature.averages$Day=='UW']
  df = subset(feature.averages,subset = Day != 'UW',invert=TRUE)
  df$diff = df[[feature.interest]]-uw_temp
  colnames(df)[1]<-'value'
  rm(p1)
  p1<-ggplot(df, aes(x = Day, y=Space)) + theme_void()+
    geom_tile(aes(fill = diff),color='black') +
    geom_text(aes(label = round(value, 1)))+
    scale_fill_gradient2(low = 'blue',mid = 'white',high = 'red',midpoint = 0)+
    labs(title = paste('UW:',as.character(round(uw_temp,1)),sep = ' '))+coord_fixed()+
    theme(axis.text = element_text(size=12),axis.title = element_text(size=12),axis.title.y = element_text(size=12,angle = 90),legend.title = element_blank())
  
  plotlist[[i]]<-p1
  
}
grid.arrange(grobs = plotlist, ncol = 6) 

#####print top contributing genes
for (i in 1:ncol(w)){
  print(i)
  print((sort(w[,i], index.return=TRUE, decreasing=TRUE)$x) %>% head(30))
}