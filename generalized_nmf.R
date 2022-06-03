library(NMF)
#calculate the most variable genes in the SCT assay
cd45.pos.sub<-FindVariableFeatures(cd45.pos.sub,selection.method = "vst",nfeatures = 3000,assay = "SCT")
plot1 = VariableFeaturePlot(cd45.pos.sub,assay = 'SCT')
plot1
#gene must be present in at least 2% of cells in subset
var.genes = VariableFeatures(cd45.pos.sub,assay = 'SCT')
derp = GetAssayData(cd45.pos.sub,assay = 'RNA',slot = "counts")
derp = derp[match(var.genes,row.names(derp)),]
pass = Matrix::rowSums(derp>0)
pass = pass/ncol(cd45.pos.sub)
pass = pass[pass>0.02]
#remove mitochondrial or ribosomal reads
###################
bad.genes = c()
bad.genes = c(bad.genes,grep(pattern = "^Rps", x = names(pass)))
bad.genes = c(bad.genes,grep(pattern = "^Rpl", x = names(pass)))
bad.genes = c(bad.genes,grep(pattern = "^mt-", x = names(pass)))
bad.genes = c(bad.genes,grep(pattern = "^Hbb-", x = names(pass)))
bad.genes = c(bad.genes,grep(pattern = "^Hba-", x = names(pass)))

pass.cleaned = pass[-bad.genes]
pass.cleaned = pass.cleaned[1:1250]
###Use Seurat's builtin scaledata function:
cd45.pos.sub = ScaleData(cd45.pos.sub,assay = 'RNA',verbose = TRUE,do.center = FALSE)
scaled = GetAssayData(cd45.pos.sub,slot = "scale.data",assay = 'RNA')
var.genes.nmf <- names(pass.cleaned)
var.norm.counts = scaled[var.genes.nmf,]
############################################################
#set range of number of factors to parameter sweep through
param.sweep = c(3:7)
#iterate through parameter sweep, multi-thread with 20cores in this instance. use the ns-nmf algorithm
for (i in 1:length(param.sweep)){
  nmf.out.ns.temp = nmf(as.matrix(var.norm.counts),rank = param.sweep[i],nrun=50,seed = 12345, method = 'nsNMF', .opt = 'vp20')
  saveRDS(nmf.out.ns.temp, file = paste("~/skeletal_seurat_scaling/nmf_ns_",as.character(param.sweep[i]),'.rds',sep = ''))
}
