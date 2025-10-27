###去批次
features <- SelectIntegrationFeatures(object.list = scRNAlist4)
for (i in 1:length(scRNAlist4)){
    scRNAlist4[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist4[[i]], pattern = "^mt-")
    scRNAlist4[[i]] <- subset(scRNAlist4[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5)
    scRNAlist4[[i]] <- NormalizeData(scRNAlist4[[i]])
    scRNAlist4[[i]] <- FindVariableFeatures(scRNAlist4[[i]],selection.method = 'vst',nfeatures = 2000) # 这里我们只需要2000的基因
    scRNAlist4[[i]] <- ScaleData(scRNAlist4[[i]], features = features, verbose = FALSE)
    scRNAlist4[[i]] <- RunPCA(scRNAlist4[[i]], features = features, verbose = FALSE)
}

# 寻找Integration anchors
# 我们使用PCA来寻找anchor
scRNA4.anchors <- FindIntegrationAnchors(object.list = scRNAlist4, dims = 1:30)
# 整合数据
scRNA4 <- IntegrateData(anchorset = scRNA4.anchors, dims = 1:30)
# 设置默认使用数据
DefaultAssay(scRNA4) <- "integrated"

# Run the standard workflow for visualization and clustering
# 标准化
scRNA4 <- ScaleData(scRNA4, verbose = FALSE)
# 使用主成分分析（PCA）执行线性降维
scRNA4 <- RunPCA(scRNA4, npcs = 30, verbose = FALSE)
#  Non-linear dimensionality reduction using UMAP (Uniform Manifold Approximation and Projection).
scRNA4  <- RunUMAP(scRNA4 , dims=1:30,reduction='pca')
# Constructs a shared nearest neighbor (SNN) graph for clustering.
scRNA4 <- FindNeighbors(scRNA4, reduction = "pca", dims = 1:30)
#  Identifies cell clusters using the Louvain algorithm (default).
scRNA4 <- FindClusters(scRNA4,resolution = 0.8)

saveRDS(scRNA4,"scRNA4.filt.anchors.rds")

# individual clusters
p1=DimPlot(scRNA4,group.by = "orig.ident")
ggsave(filename = "p_DimPlotAnchors1_scRNA4.pdf", plot = p1, height = 6, width = 6)

p2=DimPlot(scRNA4,group.by = "seurat_clusters",label=T)
ggsave(filename = "p_DimPlotAnchors2_scRNA4.pdf", plot = p2, height = 6, width = 6)

combine<-CombinePlots(list(p1,p2))
ggsave(filename = "p_DimPlotAnchors_scRNA4.pdf", plot = combine, height = 6, width = 12)

# Find marker genes
scRNA4.markers <- FindAllMarkers(scRNA4, only.pos = TRUE)
library(openxlsx)
write.xlsx(scRNA4.markers,'scRNA4.markers.xlsx')

# Cell type annotation
Adipocyte_Adipose_Tissue=c(0,1,10,15,24,29)
B_cell=c(7,19,26)
Macrophage=c(25)
Endothelial_Cell =c(3,6,22)
Fibroblast=c(4,12,16,21)
MP=c(2,20,27)
Monocyte_Cell =c(33)
Muscle_Cell=c(31,34)
NK_Cell=c(17)
myoFB=c(8,9,11,14)
Smooth_Muscle_Cell=c(13)
T_cell=c(5)
unknown=c(18,28,30,32,23)


DefaultAssay(scRNA4) <- "RNA"
# marker
genes_to_check = c("Adipoq","Pnpla2","Plin1","Cidec","Apoc1","Fabp4",
                  "Pxk","Ms4a1","Cd19","Cd74","Cd79a","Ighd",
                  "Cd83","Cd86","Ly75",
                  "Itgax","Zbtb46","Cd86","Cd83","Cd1a",
                  "Cd93","Vwf","Emcn","Egfl7","Flt1","Id3",
                  "Vim","Pdgfrb","Lum","Col6a2","Vtn","Mfap5",
                  "Cd68","Fcgr1","Lyz2","Sepp1","Naaa","Ccr2","Cd74",
                  "Tnnt3","Ttn",
                  "Map2",
                  "Nkg7","Klrf1","Klrd1","Gnly","Ncr1",
                  "Acta2","Myl9","Rgs5","Mylk","Nebl","Myh11",
                  "Trbc2","Cd3d","Cd3g","Cd3e","Il7r","Ltb",
                  "Cd34","Cd31","Cd133","Vegfr2","Vwf")


scRNA4@meta.data$seurat_clusters=factor(scRNA4@meta.data$seurat_clusters,levels =c(0,1,10,15,24,29,7,19,26,25,3,6,22,4,12,16,21,2,20,27,33,31,34,17,8,9,11,14,13,5,18,28,30,32,23))
pdf(file="p_clustersgenedotplot2.pdf",height = 12, width = 12)
DotPlot(scRNA4,group.by = 'seurat_clusters', features = unique(genes_to_check)) + coord_flip()
dev.off()

current.cluster.ids <- c(Adipocyte_Adipose_Tissue,B_cell,Endothelial_Cell,Fibroblast,Macrophage,Monocyte_Cell,MP,Muscle_Cell,myoFB,NK_Cell,Smooth_Muscle_Cell,T_cell,unknown)
new.cluster.ids <- c(rep("AP",length(Adipocyte_Adipose_Tissue)),
                     rep("B",length(B_cell)),
                     rep("EC",length(Endothelial_Cell)),
                     rep("FB",length(Fibroblast)),
                     rep("mDC",length(Macrophage)),
                     rep("Mono",length(Monocyte_Cell)),
                     rep("MP",length(MP)),
                     rep("muscle",length(Muscle_Cell)),
                     rep("myoFB",length(myoFB)),
                     rep("NK",length(NK_Cell)),
                     rep("SMC",length(Smooth_Muscle_Cell)), 
                     rep("T",length(T_cell)),
                     rep("Unknown",length(unknown))
                     )
scRNA4@meta.data$cellType <- plyr::mapvalues(x = scRNA4$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)
scRNA4@meta.data$cellType=factor(scRNA4@meta.data$cellType,levels =c("AP","B","mDC","EC","FB","MP","Mono","muscle","NK","myoFB","SMC","T","Unknown"))
saveRDS(scRNA4,"scRNA4.celltype.rds")

# individual cell type
pdf(file="p_cellType_anno_scRNA4.pdf",height = 6, width = 6)
DimPlot(scRNA4,group.by = "cellType",label=T)
dev.off()

pdf(file="p_seurat_clusters_scRNA4.pdf",height = 6, width = 6)
DimPlot(scRNA4,group.by = "seurat_clusters",label=T)
dev.off()

pdf(file="anno_p_genedotplot-cellType.pdf",height = 12, width = 6)
DotPlot(scRNA4,group.by = 'cellType', features = unique(genes_to_check)) + coord_flip() +RotatedAxis()+ theme(text = element_text(size = 16))+theme(axis.text.x = element_text(size=20))+theme(axis.text.y = element_text(size=20))
dev.off()


####各细胞类型占比####

scRNA4@meta.data$orig.ident=factor(scRNA4@meta.data$orig.ident,levels =c("Je","Se","Ji","Si"))
Ratio <- scRNA4@meta.data %>%group_by(orig.ident,cellType) %>%
  count() %>%
  group_by(orig.ident) %>%
  mutate(Freq = n/sum(n)*100)
write.csv(file='ratio.csv',Ratio)

# 不添加label
p1=ggplot(Ratio, aes(x = orig.ident, y = Freq, fill = cellType))+
  geom_col()+
  theme_classic()+
  scale_fill_manual(values = c("#FB8072", "#1965B0", "#7BAFDE", "#882E72","#B17BA6", 
                              "#FF7F00", "#FDB462", "#E7298A", "#E78AC3","#33A02C", 
                              "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D","#E6AB02"))
ggsave(filename = "p_Rationolabel.pdf", plot = p1, height = 6, width = 6)
