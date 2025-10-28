setwd('/path/publicdataana/')
all_files<- list.files('./data/')
library(Seurat)
library(DoubletFinder)
library(dplyr)
library(ggplot2)
library(sctransform)

sobj_list<- list()
for (i in 1:length(all_files)) {
  dir1<- paste0('./data/',all_files[i])
  matrix<- Read10X(dir1)
  dir2<- paste0('./CELLBENDER/',all_files[i],'_cell_barcodes.csv')
  bar<- read.csv(dir2)
  temp<- matrix[,bar[,1]]
  dd<- CreateSeuratObject(temp,project = all_files[i])
  dd[['percent.mt']]<- PercentageFeatureSet(dd,pattern = '^MT-')
  dd[['percent.ribo']]<- PercentageFeatureSet(dd,pattern = "^RP[SL]")
  dd[['log10GenesPerUMI']] <- log10(dd$nFeature_RNA) / log10(dd$nCount_RNA)
  sobj_list[[i]]<- dd
  print(i)
}
head(dd@meta.data)
for (i in 1:length(sobj_list)) {
  temp<- sobj_list[[i]]
  temp<- subset(temp,subset = nCount_RNA>=1000 &
                  nFeature_RNA>=400 &
                  percent.mt<=5 &
                  percent.ribo<=5 &
                  log10GenesPerUMI>=0.85
                  )
  sobj_list[[i]]<- temp
  print(sobj_list[[i]])
}
head(temp@meta.data)
VlnPlot(temp,features = c('nCount_RNA','nFeature_RNA' ,'percent.mt','percent.ribo','log10GenesPerUMI'))

pc.num <- 1:30
DoubletRate = 0.023 
#options(future.globals.maxSize= 850*1024^5)
for (i in 1:length(sobj_list)) {
  temp<- sobj_list[[i]]
  temp<- SCTransform(temp,vars.to.regress = 'nCount_RNA')
  temp <- RunPCA(temp, verbose = T)
  temp <- RunUMAP(temp, dims = 1:30, verbose = T)
  
  temp <- FindNeighbors(temp, dims = 1:30, verbose = FALSE)
  temp <- FindClusters(temp, verbose = FALSE)
  DimPlot(temp, label = TRUE) + NoLegend()
  
  sweep.res <- paramSweep_v3(temp, PCs = pc.num, sct = T) 
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
  
  homotypic.prop <- modelHomotypic(temp$seurat_clusters)   
  nExp_poi <- round(DoubletRate * ncol(temp))
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  
  temp <- doubletFinder_v3(temp,
                           PCs = pc.num,
                           pN = 0.25,
                           pK = pK_bcmvn, 
                           nExp = nExp_poi.adj,
                           reuse.pANN = F,
                           sct = T)
  head(temp)
  sobj_list[[i]]<- temp
  print(i)
}
head(sobj_list[[18]]@meta.data)
table(sobj_list[[18]]$DF.classifications_0.25_0.13_119)
DimPlot(sobj_list[[18]],group.by = 'DF.classifications_0.25_0.13_119')
sobj_list[[18]]<- subset(sobj_list[[18]],subset = DF.classifications_0.25_0.13_119=='Singlet')

for (i in 1:length(sobj_list)) {
  temp<- sobj_list[[i]]
  colnames(temp@meta.data)[11:12]<- c('pANN','DF.classification')
  sobj_list[[i]]<- temp
  print(table(sobj_list[[i]]$DF.classification))
}

sce <- merge(sobj_list[[1]],y=c(sobj_list[[2]],sobj_list[[3]],sobj_list[[4]],sobj_list[[5]],sobj_list[[6]],
                                sobj_list[[7]],sobj_list[[8]],sobj_list[[9]],sobj_list[[10]],sobj_list[[11]],
                                sobj_list[[12]],sobj_list[[13]],sobj_list[[14]],sobj_list[[15]],sobj_list[[16]],
                                sobj_list[[17]],sobj_list[[18]]),project = "AT",add.cell.ids =all_files)
sce
head(sce@meta.data)
table(sce$orig.ident)
sce$group<- NA
sce@meta.data[sce$orig.ident %in% c('GSM8955403','GSM8955406','GSM8955409',
                                    'GSM8955412','GSM8955415','GSM8955418'),]$group<- 'Obese'
sce@meta.data[sce$orig.ident %in% c('GSM8955404','GSM8955407','GSM8955410',
                                    'GSM8955413','GSM8955416','GSM8955419'),]$group<- 'Lean'
sce@meta.data[sce$orig.ident %in% c('GSM8955405','GSM8955408','GSM8955411',
                                    'GSM8955414','GSM8955417','GSM8955420'),]$group<- 'Weight Loss'
table(sce$group)

Idents(sce)<-sce$orig.ident
head(sce@meta.data)
table(sce$orig.ident)

sce<- SCTransform(sce)
sce <- RunPCA(sce, verbose = T)
sce <- RunUMAP(sce, dims = 1:30, verbose = T)

sce <- FindNeighbors(sce, dims = 1:30, verbose = T)
sce <- FindClusters(sce, verbose = T)

DimPlot(sce,group.by = 'seurat_clusters',reduction = 'umap',label = T) |
  DimPlot(sce,group.by = 'orig.ident',reduction = 'umap')

save(sobj_list,file = 'sobj_list.RData')
save(sce,file = 'sce.RData')

####
sobj.anchors <- FindIntegrationAnchors(object.list = sobj_list, dims = 1:30)
sobj <- IntegrateData(anchorset = sobj.anchors, dims = 1:30)
DefaultAssay(sobj) <- "integrated"
sobj<- SCTransform(sobj)
sobj <- RunPCA(sobj, npcs = 45, verbose = FALSE)
sobj <- FindNeighbors(sobj, reduction = "pca", dims = 1:45)
sobj <- FindClusters(sobj,resolution = 0.15)
sobj <- RunUMAP(sobj, reduction = "pca", dims = 1:45)

DimPlot(sobj,group.by = 'seurat_clusters',reduction = 'umap',label = T,raster=FALSE) 
DimPlot(sobj,group.by = 'orig.ident',reduction = 'umap',raster=FALSE)

save(sobj,file = 'sobj-completed.RData')

markers_list<- c(
  'PTPRC',#IMMUNE
  'PPARG','ADIPOQ','GPAM','PLIN4',#Adipocytes #DE9B13
  'IGHM','BANK1','SELL',#B cells #56ABDB
  'VWF','BTNL9','CDH5',#Endothelial #10986F
  'MFAP5','CD55','ITGA11','DCN','LUM','PDGFRA',#ASC #E6DC48
  #'DCN','LUM','PDGFRA',#APC #0A6CA7
  'PROX1','MMRN1','CCL21',#Lymphatic #CC5D17
  'MRC1','CD163','MSR1',#Macrophages #C2749F
  'CIITA','FLT3','FCN1',#Monos/DCs #666666
  'CPA3','KIT','IL1RL1',#Mast #A87522
  'POSTN','RGS5','KCNJ8',#Mural #2488C5
  'NKG7','KLRF1','NCAM1',#NK cells #047353
  'CD3E','CD4','LEF1',#CD4+ T cells #CDC021
  'CD8A','CD8B','DTHD1',#CD8+ T cells #02537F
  'AFF3','RUNX2','ATP8B4' #ILCs #994723
)
DotPlot(sobj,features = markers_list)+coord_flip()

sobj$celltype<- as.character(sobj$seurat_clusters)
table(sobj$celltype)

sobj@meta.data[sobj$seurat_clusters %in% c(1,3,6,15),]$celltype<- 'Adipocytes'
sobj@meta.data[sobj$seurat_clusters %in% c(10),]$celltype<- 'B cells'
sobj@meta.data[sobj$seurat_clusters %in% c(4),]$celltype<- 'Endothelial'
sobj@meta.data[sobj$seurat_clusters %in% c(0),]$celltype<- 'ASC/APC'
sobj@meta.data[sobj$seurat_clusters %in% c(11),]$celltype<- 'ASC/APC'
sobj@meta.data[sobj$seurat_clusters %in% c(12),]$celltype<- 'Lymphatic'
sobj@meta.data[sobj$seurat_clusters %in% c(2,5,13),]$celltype<- 'Macrophages'
sobj@meta.data[sobj$seurat_clusters %in% c(14),]$celltype<- 'Monos/DCs'
sobj@meta.data[sobj$seurat_clusters %in% c(9),]$celltype<- 'Mast'
sobj@meta.data[sobj$seurat_clusters %in% c(8),]$celltype<- 'Mural'

sobj@meta.data[sobj$seurat_clusters %in% c(7),]$celltype<- 'T/NK'
sobj@meta.data[sobj$seurat_clusters %in% c(16),]$celltype<- 'ILCs Kit'

DimPlot(sobj,group.by = 'celltype',label = T,raster = F)
pdf(file="PubData_celltype_anno_sobj.pdf",height = 6, width = 6)
DimPlot(sobj,group.by = "celltype",label=T,raster = F)
dev.off()
DotPlot(sobj,features = markers_list,group.by = 'celltype')+coord_flip()

save(sobj,file = 'sobj-completed.RData')

sobj$group<- ''
sobj@meta.data[sobj$orig.ident %in% c('GSM8955403','GSM8955406','GSM8955409','GSM8955412','GSM8955415','GSM8955418'),]$group<- 'Obese'
sobj@meta.data[sobj$orig.ident %in% c('GSM8955404','GSM8955407','GSM8955410','GSM8955413','GSM8955416','GSM8955419'),]$group<- 'Lean'
sobj@meta.data[sobj$orig.ident %in% c('GSM8955405','GSM8955408','GSM8955411','GSM8955414','GSM8955417','GSM8955420'),]$group<- 'Weight loss'

table(sobj$group)
summary(table(sobj$orig.ident))

library(progeny)
library(readr)
library(pheatmap)
library(tibble)
library(GSVA)
sobjOW=subset(sobj,subset = group=='Obese'|group=='Weight loss')
sobjOWEC=subset(x=sobjOW,subset=(celltype == "Endothelial"))

library(clusterProfiler)
genelist<- list.files('./gmtps_hs/')
fgsea_sets<- list()
for (i in 1:length(genelist)) {
  dir<- paste0('./gmtps_hs/',genelist[i])
  temp<- read.gmt(dir)
  fgsea_sets[[names(table(temp$term))]]<- temp$gene
}
sobjOWEC <- AddModuleScore(
  object = sobjOWEC,
  features = fgsea_sets[5],
  name = "VEGF"
)

head(sobjOWEC@meta.data)
colnames(sobjOWEC@meta.data)[length(colnames(sobjOWEC@meta.data))]='REACTOME_VEGF_LIGAND_RECEPTOR_INTERACTIONS'
library(ggpubr)
ggviolin(sobjOWEC@meta.data,
         x="group",y="REACTOME_VEGF_LIGAND_RECEPTOR_INTERACTIONS",
         width=0.8,color="black",
         fill="group",
         xlab=F,
         add='mean_sd',
         bxp.errorbar=T,
         bxp.errorbar.width=0.05,
         size=0.5,
         palette="npg",
         legend="right")
ggsave(filename="genesetsScoreviolin.pdf",width=4,height=3)


