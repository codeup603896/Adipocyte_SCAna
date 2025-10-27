# brown AP

scRNA4SeJe=subset(x = scRNA4, subset = (orig.ident == "Se"|orig.ident == "Je"))

scRNA4SeJeAP=subset(x = scRNA4SeJe, subset = (cellType == "AP"))
# AP亚群分类
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


# 棕色脂肪拟时序分析
library(monocle)
library(RColorBrewer)

# AP+myoFB细胞分化轨迹
scRNA4SeJeAPmyoFB=subset(x = scRNA4SeJe, subset = (cellType == "AP"|cellType == "myoFB"))

scRNA4SeJeAP@meta.data$cellType2<- paste0(scRNA4SeJeAP@meta.data$cellType,scRNA4SeJeAP@meta.data$seurat_clusters)
library(dplyr)
 
# 将df1转换为tibble，以便使用left_join，并设置行名为一个临时列
scRNA4SeJeAP_tibble <- tibble::rownames_to_column(as_tibble(scRNA4SeJeAP@meta.data), var = "rowname")
scRNA4SeJeAP_tibble$rowname=rownames(scRNA4SeJeAP@meta.data)
scRNA4SeJeAPmyoFB_tibble <- tibble::rownames_to_column(as_tibble(scRNA4SeJeAPmyoFB@meta.data), var = "rowname")
scRNA4SeJeAPmyoFB_tibble$rowname=rownames(scRNA4SeJeAPmyoFB@meta.data)
# 使用left_join合并数据框，然后根据结果替换df2中的celltype2列
merged_df <- left_join(scRNA4SeJeAPmyoFB_tibble, scRNA4SeJeAP_tibble, by = "rowname")
 
# 检查是否有NA值，如果有，则用原始的celltype2值填充
merged_df$cellType2 <- as.character(merged_df$cellType2)
merged_df$cellType.x<- as.character(merged_df$cellType.x)
# 使用ifelse进行替换
merged_df$cellType2 <- ifelse(is.na(merged_df$cellType2), merged_df$cellType.x, merged_df$cellType2)

scRNA4SeJeAPmyoFB@meta.data$cellType2=merged_df$cellType2

expr_matrix <- as.data.frame(scRNA4SeJeAPmyoFB@assays$RNA@counts)
# 细胞注释矩阵（列为细胞名）
sample_sheet <- scRNA4SeJeAPmyoFB@meta.data
# 基因注释矩阵（行为基因名）
gene_annotation <- data.frame(
  gene_short_name=rownames(scRNA4SeJeAPmyoFB@assays$RNA),
  row.names = rownames(scRNA4SeJeAPmyoFB@assays$RNA)
)

pd <- new("AnnotatedDataFrame", data = sample_sheet)
fd <- new("AnnotatedDataFrame", data = gene_annotation)

#创建monocle对象
cd <- newCellDataSet(as(as.matrix(expr_matrix),"sparseMatrix"), 
                     phenoData = pd, 
                     featureData = fd,
                     lowerDetectionLimit = 0.5,
                     expressionFamily = negbinomial.size())

###计算数据的size factors和dispersions
cd <- estimateSizeFactors(cd)
cd <- estimateDispersions(cd)
# 过滤低于1%细胞中检出的基因，最低表达阈值为0.5
cd <- detectGenes(cd, min_expr = 0.5)
expressed_genes <- row.names(subset(fData(cd), num_cells_expressed > nrow(sample_sheet) * 0.01))
# 查看筛选后的基因个数（用于后面基因的筛选）
length(expressed_genes)

###构建轨迹
disp_table <- dispersionTable(cd)
# 高度离散基因的筛选标准，可根据数据情况设置mean_expression的值
ordering_genes <- as.character(subset(disp_table,mean_expression >= 0.3 & dispersion_empirical >= dispersion_fit)$gene_id)
cd <- setOrderingFilter(cd, ordering_genes)
# plot_ordering_genes(cd)
# 根据数据的某种分类进行差异分析
diff_test_res <- differentialGeneTest(cd[expressed_genes,],fullModelFormulaStr = "~cellType2")
# 筛选差异基因（q < 1e-5并且属于之前计算的expressed_genes列表中）
ordering_genes <- row.names (subset(diff_test_res, qval < 1e-5))
ordering_genes <- intersect(ordering_genes, expressed_genes)
cd <- setOrderingFilter(cd, ordering_genes)
# plot_ordering_genes(cd)
# 选取的基因数目为每个细胞的维度，基于默认的'DDRTree'方法进行数据降维
cd <- reduceDimension(cd, max_components = 2, method = 'DDRTree')
# 对细胞进行排序，由于排序无法区分起点和终点，若分析所得时序与实际相反，根据“reverse”参数进行调整，默认reverse=F
cd <- orderCells(cd, reverse = F) 
# 排序好的细胞可以进行可视化，可标注细胞的各注释信息)
# 查看细胞注释信息
head(cd@phenoData@data)
saveRDS(cd,file='SeJem2_cd_apmyofb.rds')

# 设置颜色
color1 <- c(brewer.pal(8, "Set1"))
getPalette <- colorRampPalette(brewer.pal(6, "Set1"))

# 拆分
pData(cd)$sample=scRNA4SeJeAPmyoFB$orig.ident
pData(cd)$sample=factor(pData(cd)$sample, levels = c('Je','Se'))
pData(cd)$cellType2 <- factor(pData(cd)$cellType2, levels = c('AP4','AP2','AP1','AP7','AP0','AP5','AP3','AP6','myoFB'))


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
