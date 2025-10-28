# brown AP

scRNA4SeJe=subset(x = scRNA4, subset = (orig.ident == "Se"|orig.ident == "Je"))
scRNA4SeJeAP=subset(x = scRNA4SeJe, subset = (cellType == "AP"))

DefaultAssay(scRNA4SeJeAP) <- "RNA"
scRNA4SeJeAP <- NormalizeData(scRNA4SeJeAP, normalization.method = "LogNormalize", scale.factor = 10000)
scRNA4SeJeAP <- FindVariableFeatures(scRNA4SeJeAP,nfeatures = 2000)
scRNA4SeJeAP <- ScaleData(scRNA4SeJeAP,vars.to.regress = 'percent.mt')
# scRNA4SeJeAP <- ScaleData(scRNA4SeJeAP, verbose = FALSE)
scRNA4SeJeAP <- RunPCA(scRNA4SeJeAP)
scRNA4SeJeAP <- FindNeighbors(scRNA4SeJeAP, reduction = "pca", dims = 1:30)
scRNA4SeJeAP <- FindClusters(scRNA4SeJeAP, resolution = 0.3)
scRNA4SeJeAP <- RunUMAP(scRNA4SeJeAP, reduction = "pca", dims = 1:30) 
library(harmony)
scRNA4SeJeAP <- RunHarmony(scRNA4SeJeAP, group.by.vars = "orig.ident")
scRNA4SeJeAP  <- FindNeighbors(scRNA4SeJeAP , dims = 1:30 , reduction = "harmony")
scRNA4SeJeAP  <- FindClusters(scRNA4SeJeAP, save.snn=T , resolution = 0.3)
scRNA4SeJeAP  <- RunUMAP(scRNA4SeJeAP , dims=1:30,reduction='harmony')
table(scRNA4SeJeAP@meta.data$seurat_clusters)

pdf(file=sprintf("scRNA4SeJeAP_Nucb2.dotplot.pdf"),height = 6, width = 6)
DotPlot(scRNA4SeJeAP,group.by = 'seurat_clusters', features = c('Nucb2')) + coord_flip()
dev.off()
saveRDS(scRNA4SeJeAP,'scRNA4SeJeAP_harmony.rds')


library(monocle)
library(RColorBrewer)
scRNA4SeJeAPFBO=subset(x = scRNA4SeJe, subset = (cellType == "AP"|cellType == "FBO"))
scRNA4SeJeAP@meta.data$cellType2<- paste0(scRNA4SeJeAP@meta.data$cellType,scRNA4SeJeAP@meta.data$seurat_clusters)
library(dplyr)

scRNA4SeJeAP_tibble <- tibble::rownames_to_column(as_tibble(scRNA4SeJeAP@meta.data), var = "rowname")
scRNA4SeJeAP_tibble$rowname=rownames(scRNA4SeJeAP@meta.data)
scRNA4SeJeAPFBO_tibble <- tibble::rownames_to_column(as_tibble(scRNA4SeJeAPFBO@meta.data), var = "rowname")
scRNA4SeJeAPFBO_tibble$rowname=rownames(scRNA4SeJeAPFBO@meta.data)
merged_df <- left_join(scRNA4SeJeAPFBO_tibble, scRNA4SeJeAP_tibble, by = "rowname")
merged_df$cellType2 <- as.character(merged_df$cellType2)
merged_df$cellType.x<- as.character(merged_df$cellType.x)
merged_df$cellType2 <- ifelse(is.na(merged_df$cellType2), merged_df$cellType.x, merged_df$cellType2)
scRNA4SeJeAPFBO@meta.data$cellType2=merged_df$cellType2

expr_matrix <- as.data.frame(scRNA4SeJeAPFBO@assays$RNA@counts)
sample_sheet <- scRNA4SeJeAPFBO@meta.data
gene_annotation <- data.frame(
  gene_short_name=rownames(scRNA4SeJeAPFBO@assays$RNA),
  row.names = rownames(scRNA4SeJeAPFBO@assays$RNA)
)
pd <- new("AnnotatedDataFrame", data = sample_sheet)
fd <- new("AnnotatedDataFrame", data = gene_annotation)
cd <- newCellDataSet(as(as.matrix(expr_matrix),"sparseMatrix"), 
                     phenoData = pd, 
                     featureData = fd,
                     lowerDetectionLimit = 0.5,
                     expressionFamily = negbinomial.size())

cd <- estimateSizeFactors(cd)
cd <- estimateDispersions(cd)
cd <- detectGenes(cd, min_expr = 0.5)
expressed_genes <- row.names(subset(fData(cd), num_cells_expressed > nrow(sample_sheet) * 0.01))

disp_table <- dispersionTable(cd)
ordering_genes <- as.character(subset(disp_table,mean_expression >= 0.3 & dispersion_empirical >= dispersion_fit)$gene_id)
cd <- setOrderingFilter(cd, ordering_genes)
# plot_ordering_genes(cd)

diff_test_res <- differentialGeneTest(cd[expressed_genes,],fullModelFormulaStr = "~cellType2")
ordering_genes <- row.names (subset(diff_test_res, qval < 1e-5))
ordering_genes <- intersect(ordering_genes, expressed_genes)
cd <- setOrderingFilter(cd, ordering_genes)
# plot_ordering_genes(cd)

cd <- reduceDimension(cd, max_components = 2, method = 'DDRTree')
cd <- orderCells(cd, reverse = F) 
saveRDS(cd,file='SeJem2_cd_apFBO.rds')

color1 <- c(brewer.pal(8, "Set1"))
getPalette <- colorRampPalette(brewer.pal(6, "Set1"))

pData(cd)$sample=scRNA4SeJeAPFBO$orig.ident
pData(cd)$sample=factor(pData(cd)$sample, levels = c('Je','Se'))
pData(cd)$cellType2 <- factor(pData(cd)$cellType2, levels = c('AP4','AP2','AP1','AP7','AP0','AP5','AP3','AP6','FBO'))

pdf(file="sejecd_cellType.trajectory.pdf",width=12,height=6)
plot_cell_trajectory(cd, color_by = "cellType", cell_size = 0.5,cell_link_size = 0.5) + facet_wrap(~sample, nrow = 1, scales = "free")+ scale_color_manual(values=c('red','blue'))+ggtitle("cellType") + theme(plot.title = element_text(hjust = 0.5),legend.position = "right")
dev.off()

pdf(file="sejecd_sample.trajectory.pdf",width=12,height=6)
plot_cell_trajectory(cd, color_by = "sample", cell_size = 0.5,cell_link_size = 0.5) + facet_wrap(~cellType2, nrow = 3, scales = "free")+ggtitle("sample") + theme(plot.title = element_text(hjust = 0.5),legend.position = "right")
dev.off()

BEAM_res <- BEAM(cd, branch_point = 3, cores = 50, progenitor_method = "duplicate")
write.csv(file="SeJe.BEAM_res.csv",BEAM_res,rowname)
# BEAM_res=read.csv(file="SeJe.BEAM_res.csv")
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
pdf(file="sejem2_genes_branched_heatmap50.pdf",width=4,height=6)
plot_genes_branched_heatmap(cd[c(row.names(subset(BEAM_res[1:50,],qval < 1e-4))),],
                            branch_point = 3,
                            num_clusters = 4,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T)
dev.off()

