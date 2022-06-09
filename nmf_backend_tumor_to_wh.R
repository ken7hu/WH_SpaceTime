# Script for comparing factor similarity between wound healing and tumor settings
nmf.out.wh = readRDS("~/Krummel lab/20201130_WH011/WH_REDO/nmf/cd45pos/second_try/mm_outs_seurat_scaling/nmf_ns_24.rds")
nmf.out.tumor = readRDS("~/Krummel lab/20201130_WH011/WH_REDO/nmf/humu/nmf_seurat_scaling_2nd/nmf_out_25.rds")
tumor_w = basis(nmf.out.tumor)
wh_w = basis(nmf.out.wh)
mutual_genes = intersect(rownames(wh_w),rownames(tumor_w))
####specify factors you want plotted against each other, for the scatterplots Figure 6C-F
tumor_factor = 12
wh_factor =1
weights = data.frame(mutual_genes)
temp = wh_w[,wh_factor]
weights$wh_weights = temp[match(mutual_genes,names(temp))]
temp = tumor_w[,tumor_factor]
weights$tumor_weights = temp[match(mutual_genes,names(temp))]
###for labeling the points
n.label = 20
weights$label = rep('',nrow(weights))
lst <- sort(weights$wh_weights, index.return=TRUE, decreasing=TRUE)
cutoff_1 = lst$x[n.label]
labels = lst$ix %>% head(n.label)
lst <- sort(weights$tumor_weights, index.return=TRUE, decreasing=TRUE)
labels = union(labels,lst$ix %>% head(n.label))

cutoff_2 = lst$x[n.label]
weights$label[labels] = as.character(weights$mutual_genes[labels])


g6 = ggplot(data = weights,aes(x=`wh_weights`,y=`tumor_weights`,label = label))+geom_point()+
  geom_text_repel()+geom_abline(slope=1)+theme_classic()+geom_hline(yintercept = cutoff_2,linetype = 'dashed')+
  geom_vline(xintercept = cutoff_1,linetype = 'dashed')+ylab(paste('Gene Weight Tumor Factor ',as.character(tumor_factor),sep=''))+xlab(paste('Gene Weight WH Factor ',as.character(wh_factor),sep=''))+
  theme(axis.title = element_text(size=15),axis.text = element_text(size=14))
g6

g1+g2+g3+g4+g5+g6+plot_layout(ncol = 2)

#### heatmap with jaccard dist. between n tumor and m wound healing factors
corr = matrix(0,ncol(wh_w),ncol(tumor_w))
alpha = matrix(0,ncol(wh_w),ncol(tumor_w))
row.names(corr) = paste('factor-',c(1:ncol(wh_w)),sep='')
colnames(corr) =  paste('factor-',c(1:ncol(tumor_w)),sep='')
n.label = 20
for(i in 1:ncol(wh_w)){
  temp_wh = wh_w[,i]
  temp_wh = temp_wh[match(mutual_genes,names(temp_wh))]
  for(j in 1:ncol(tumor_w)){
    temp_tumor = tumor_w[,j]
    temp_tumor = temp_tumor[match(mutual_genes,names(temp_tumor))]
    
    lst <- sort(temp_wh, index.return=FALSE, decreasing=TRUE)
    labels = lst %>% head(n.label) %>% names()
    lst <- sort(temp_tumor, index.return=FALSE, decreasing=TRUE)
    labels2 = lst %>% head(n.label) %>% names()
    jaccard = length(intersect(labels,labels2))/length(union(labels,labels2))
    corr[i,j] = jaccard
  }
}
library(RColorBrewer)
library(viridis)
library(ComplexHeatmap)
ht2 = Heatmap(corr, name = "ht2",col = inferno(20),cluster_rows = FALSE,cluster_columns = FALSE)
hts
