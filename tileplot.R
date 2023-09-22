#TILEPLOTTING
####TILEPLOT for metadata#####
#Note most of these functions were designed to work with Seurat v3 which handles the output of AverageExpression differently from v4
TilePlot_meta <- function(obj, feature){
  require(plotly)
  if (feature %in% colnames(obj@meta.data)){
    ids = (table(obj$PunchID))
    obj$Punch_ID <- factor(obj$PunchID, levels = names(ids))
    Idents(obj)<-obj$Punch_ID
    feature.averages = aggregate(obj[[feature]], list(obj$Punch_ID), mean)
    blah = as.matrix(t(feature.averages$NMF))
    feature.averages$Day = sub("\\_.*", "", feature.averages$Group.1)
    feature.averages$Space = sub(".*\\_","",feature.averages$Group.1)
    uw_temp = feature.averages[[feature]][feature.averages$Day=='UW']
    df = subset(feature.averages,subset = Day != 'UW',invert=TRUE)
    df$diff = df[[feature]]-uw_temp
    df$val = df[[feature]]
    ggplot(df, aes(x = Day, y=Space)) +
      geom_tile(aes(fill = diff),color='black') +
      geom_text(aes(label = round(val, 3)))+
      scale_fill_gradient2(low = 'blue',mid = 'white',high = 'red',midpoint = 0)+
      labs(title = paste(feature,'UW :',as.character(round(uw_temp,3)),sep = ' '))
  }
  else{
    print("Missing")
  }
}
#################
#THIS FUNCTION PLOTS TILES FOR A GENE
#Note this function is designed to work with Seurat v3
TilePlot_gene <- function(obj, feature,assay,slot,scale.to.uw){
  ref.data = GetAssayData(obj,assay=assay,slot = "counts")
  if (feature %in% rownames(ref.data)){
    ids = (table(obj$PunchID))
    obj$Punch_ID <- factor(obj$PunchID, levels = names(ids))
    Idents(obj)<-obj$Punch_ID
    feature.averages = AverageExpression(obj,assays = assay,slot = slot,features = feature,verbose = 'FALSE')
    #feature.averages = aggregate(obj[[feature]], list(obj$Punch_ID), mean)
    feature.averages = as.data.frame(t(feature.averages[[1]]))
    feature.averages$Day = sub("\\_.*", "", row.names(feature.averages))
    feature.averages$Space = sub(".*\\_","",row.names(feature.averages))
    uw_temp = feature.averages[[feature]][feature.averages$Day=='UW']
    df = subset(feature.averages,subset = Day != 'UW',invert=TRUE)
    if(scale.to.uw){
      df$diff = df[[feature]]-uw_temp
    }
    else{
      df$diff = df[[feature]]
    }
    df$val = df[[feature]]
    g <- ggplot(df, aes(x = Day, y=Space)) +
      geom_tile(aes(fill = diff),color='black') +
      geom_text(aes(label = round(val, 1)),size=6)+
      scale_fill_gradient2(low = 'blue',mid = 'white',high = 'red',midpoint = 0,na.value = 'grey')+theme_void()+
      labs(title = paste('UW :',as.character(round(uw_temp,1)),sep = ' '))+coord_fixed()+
      theme(axis.text = element_text(size=12),axis.title = element_text(size=12),axis.title.y = element_text(size=12,angle = 90),legend.title = element_blank(),
            legend.text = element_text(size=12),title = element_text(size = 16))
  }
  else{
    print("METADATA IS  MISSING")
  }
  return(g)
}

#################
#THIS FUNCTION PLOTS TILES FOR A GENE
#Note this code is designed to work with Seurat v4
TilePlot_gene <- function(obj, feature,assay,slot,scale.to.uw){
  ref.data = GetAssayData(obj,assay=assay,slot = "counts")
  if (feature %in% rownames(ref.data)){
    ids = (table(obj$PunchID))
    obj$Punch_ID <- factor(obj$PunchID, levels = names(ids))
    Idents(obj)<-obj$Punch_ID
    feature.averages = AverageExpression(obj,assays = assay,slot = slot,features = feature,verbose = 'FALSE')
    #feature.averages = aggregate(obj[[feature]], list(obj$Punch_ID), mean)
    feature.averages = as.data.frame(t(feature.averages[[1]]))
    feature.averages$Day = sub("\\_.*", "", row.names(feature.averages))
    feature.averages$Space = sub(".*\\_","",row.names(feature.averages))
    uw_temp = feature.averages$V1[feature.averages$Day=='UW']
    df = subset(feature.averages,subset = Day != 'UW',invert=TRUE)
    if(scale.to.uw){
      df$diff = df$V1-uw_temp
    }
    else{
      df$diff = df$V1
    }
    df$val = df$V1
    g <- ggplot(df, aes(x = Day, y=Space)) +
      geom_tile(aes(fill = diff),color='black') +
      geom_text(aes(label = round(val, 1)),size=6)+
      scale_fill_gradient2(low = 'blue',mid = 'white',high = 'red',midpoint = 0,na.value = 'grey')+theme_void()+
      labs(title = paste('UW :',as.character(round(uw_temp,1)),sep = ' '))+coord_fixed()+
      theme(axis.text = element_text(size=12),axis.title = element_text(size=12),axis.title.y = element_text(size=12,angle = 90),legend.title = element_blank(),
            legend.text = element_text(size=12),title = element_text(size = 16))
  }
  else{
    print("METADATA IS  MISSING")
  }
  return(g)
}
#########################

TilePlot_gene_diff <- function(obj, feature,feature2,assay,slot,scale.to.uw){
  ref.data = GetAssayData(obj,assay=assay,slot = "counts")
  if (feature %in% rownames(ref.data)){
    ids = (table(obj$PunchID))
    obj$Punch_ID <- factor(obj$PunchID, levels = names(ids))
    Idents(obj)<-obj$Punch_ID
    feature.averages = AverageExpression(obj,assays = assay,slot = slot,features = feature,verbose = 'FALSE')
    #feature.averages = aggregate(obj[[feature]], list(obj$Punch_ID), mean)
    feature.averages = as.data.frame(t(feature.averages[[1]]))
    feature.averages$Day = sub("\\_.*", "", row.names(feature.averages))
    feature.averages$Space = sub(".*\\_","",row.names(feature.averages))
    uw_temp = feature.averages[[feature]][feature.averages$Day=='UW']
    df = subset(feature.averages,subset = Day != 'UW',invert=TRUE)
    if(scale.to.uw){
      df$diff = df[[feature]]-uw_temp
    }
    else{
      df$diff = df[[feature]]
    }
    df$val = df[[feature]]
    g <- ggplot(df, aes(x = Day, y=Space)) +
      geom_tile(aes(fill = diff),color='black') +
      geom_text(aes(label = round(val, 3)))+
      scale_fill_gradient2(low = 'blue',mid = 'white',high = 'red',midpoint = 0,na.value = 'grey')+
      labs(title = paste(feature,'UW :',as.character(round(uw_temp,3)),sep = ' '))+coord_fixed()
  }
  else{
    print("METADATA IS  MISSING")
  }
  return(g)
}



TilePlot_gene_2 <- function(obj, feature,assay,rm){
  ref.data = GetAssayData(obj,assay=assay,slot = "counts")
  if (feature %in% rownames(ref.data)){
    ids = (table(obj$PunchID))
    obj$Punch_ID <- factor(obj$PunchID, levels = names(ids))
    Idents(obj)<-obj$Punch_ID
    feature.averages = AverageExpression(obj,assays = assay,slot = 'counts',features = feature,verbose = 'FALSE')
    #feature.averages = aggregate(obj[[feature]], list(obj$Punch_ID), mean)
    feature.averages = as.data.frame(t(feature.averages[[1]]))
    feature.averages$Day = sub("\\_.*", "", row.names(feature.averages))
    feature.averages$Space = sub(".*\\_","",row.names(feature.averages))
    uw_temp = feature.averages[[feature]][feature.averages$Day=='UW']
    df = subset(feature.averages,subset = Day != 'UW',invert=TRUE)
    df$diff = df[[feature]]-uw_temp
    df$val = df[[feature]]
    df$val[rm]<-NA
    ggplot(df, aes(x = Day, y=Space)) +
      geom_tile(aes(fill = val),color='black') +
      geom_text(aes(label = round(val, 3)))+
      scale_fill_gradient(low = 'blue',high = 'red',na.value = 'grey')+
      labs(title = paste(feature,'UW :',as.character(round(uw_temp,3)),sep = ' '))+coord_fixed()
  }
  else{
    print("METADATA IS  MISSING")
  }
}
######################################3
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
    geom_text(aes(label = round(percents, 3)))+
    scale_fill_gradient2(low = 'blue',mid = 'white',high = 'red',midpoint = 0)+
    labs(title = paste(unique(blah$type)[i],'UW % :',as.character(round(uw_temp,3)),sep = ' '))
  plotlist[[i]]<-p1
}
grid.arrange(grobs = plotlist, ncol = 5) 



#############FOR FREQUENCIES#############
ln$orig.ident = Idents(ln)
y = table(ln@meta.data[,c('orig.ident','PunchID')])
y.freqs = apply(y,1,"/",colSums(y))
y.freqs = as.data.frame(y.freqs)
blah = as.matrix(y.freqs$Bplasma)
z = blah[1:16,1]
z = as.numeric(matrix(z,nrow=4))
dim(z) <- c(4,4)
uw = matrix( rep( 0, len=16), nrow = 4)
uw = uw + blah[17,1]
fig <- plot_ly(showscale = FALSE)
fig <- fig %>% add_surface(z = ~z)
fig <- fig %>% add_surface(z = uw, opacity = 0.7)
fig

#####start here?
ln$orig.ident = ln$annotation_1
ln$orig.ident = Idents(ln)

y = table(ln@meta.data[,c('orig.ident','PunchID')])
y.freqs = apply(y,1,"/",colSums(y))
y.freqs = as.data.frame(y.freqs)
library(reshape2)
y.freqs$Punch_ID = row.names(y.freqs)
blah = melt(y.freqs,id.vars = 'Punch_ID')
blah$Day = sub('_.*', '', blah$Punch_ID)
blah$Space = sub('.*_', '', blah$Punch_ID)
library(gridExtra)
plotlist = list()
for (i in 1:length(unique(blah$variable))){
  
  df = subset(blah,subset = variable == unique(blah$variable)[i])
  uw_temp = df$value[df$Day=='UW']
  df = subset(df,subset = Day != 'UW',invert=TRUE)
  df$percents_diff = df$value-uw_temp
  rm(p1)
  p1<-ggplot(df, aes(x = Day, y=Space)) +
    geom_tile(aes(fill = percents_diff),color='black') +
    geom_text(aes(label = round(value, 3)))+
    scale_fill_gradient2(low = 'blue',mid = 'white',high = 'red',midpoint = 0)+
    labs(title = paste(unique(blah$variable)[i],'UW % :',as.character(round(uw_temp,3)),sep = ' '))+coord_fixed()+
    theme_void()+theme(axis.text = element_text(size = 16),legend.title = element_blank())
  plotlist[[i]]<-p1
}
grid.arrange(grobs = plotlist, ncol = 5) 
