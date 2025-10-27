# AP6 v.s. APNo6差异基因分析

scRNA4SeJeAP@meta.data$cellType2<- 'APNo6'
Idents(scRNA4SeJeAP)='seurat_clusters'
filtered_cellsAP6 <- WhichCells(scRNA4SeJeAP,idents='6')
Idents(scRNA4SeJeAP)='cellType2'
Idents(scRNA4SeJeAP, cells = filtered_cellsAP6) <- 'AP6'
table(Idents(scRNA4SeJeAP))
scRNA4SeJeAP@meta.data$cellType2=Idents(scRNA4SeJeAP)

library(limma)
library(edgeR)
data=scRNA4SeJeAP@assays$RNA@data
dd_CM<- as.data.frame(as.matrix(data))

group2<-  factor(scRNA4SeJeAP@meta.data$cellType2,levels = c("APNo6", "AP6"))#顺序靠后的Se vs 顺序靠前的Je;#对照组在前，处理组在后
table(group2)
dgelist <- DGEList(counts =dd_CM, group = group2)
#（2）过滤 low count 数据，例如 CPM 标准化（推荐）
keep <- rowSums(cpm(dgelist) > 1 ) >= 2
dgelist <- dgelist[keep, , keep.lib.sizes = FALSE]
#（3）标准化，以 TMM 标准化为例
dgelist_norm <- calcNormFactors(dgelist, method = 'TMM')
#差异表达基因分析
design <- model.matrix(~group2)
#（1）估算基因表达值的离散度
dge <- estimateDisp(dgelist_norm, design, robust = TRUE)
#（2）模型拟合，edgeR 提供了多种拟合算法
#负二项广义对数线性模型
fit <- glmFit(dge, design, robust = TRUE)
lrt <- topTags(glmLRT(fit), n = nrow(dgelist$counts))
sample1='AP6'
sample2='APNo6'

logFC <- lrt$table$logFC
lrt$table$SE <- abs(lrt$table$logFC) / sqrt(lrt$table$LR)
# 计算 95% CI（正态近似）
lrt$table$CI.L <- logFC - 1.96 * SE
lrt$table$CI.R <- logFC + 1.96 * SE
# 合并结果
results <- data.frame(logFC, CI.L, CI.R, p.value = lrt$table$PValue)

write.table(lrt, sprintf('%svs%s.txt',sample1,sample2), sep = '\t', col.names = NA, quote = FALSE)
library(ggplot2)

d1<- read.table(sprintf('%svs%s.txt',sample1,sample2))
colnames(d1)
d1$gene<- rownames(d1)

# 设置p_value和logFC的阈值
cut_off_FDR = 0.05  #统计显著性
cut_off_logFC = 0.5           #差异倍数值
d1$change<- 'Stable'
# 根据阈值参数，上调基因设置为‘up’，下调基因设置为‘Down’，无差异设置为‘Stable’，并保存到change列中
d1$change = ifelse(d1$FDR< cut_off_FDR & abs(d1$logFC) > cut_off_logFC, 
                        ifelse(d1$logFC> cut_off_logFC ,'Up','Down'),
                        'Stable')
head(d1)
table(d1$change)
# 火山图
p <- ggplot(
  # 数据、映射、颜色
  d1, aes(x = logFC, y = -log10(FDR), colour=change)) +
  geom_point(alpha=0.4, size=2) +
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
  # 辅助线
  geom_vline(xintercept=c(-0.5,0.5),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_FDR),lty=4,col="black",lwd=0.8) +
  # 坐标轴
  labs(x="log2(fold change)",
       y="-log10 (FDR)")+
  scale_x_continuous(limits = c(-5, 5))+
  scale_y_continuous(limits = c(0,30))+
  theme_bw()+
  # 图例
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())+
  ggtitle('AP6 v.s. APNo6')
# p
ggsave(filename = "diffgene2nolabel.pdf", height = 6, width = 8, plot = p)


# 获取不同类型细胞所有基因的表达平均值及表达比例
scRNA4SeJeAPmyoFBEC=subset(x = scRNA4SeJe, subset = (cellType == "AP"|cellType == "myoFB"|cellType == "EC"))

library(dplyr)
 
# 将df1转换为tibble，以便使用left_join，并设置行名为一个临时列
scRNA4SeJeAP_tibble <- tibble::rownames_to_column(as_tibble(scRNA4SeJeAP@meta.data), var = "rowname")
scRNA4SeJeAP_tibble$rowname=rownames(scRNA4SeJeAP@meta.data)
scRNA4SeJeAPmyoFBEC_tibble <- tibble::rownames_to_column(as_tibble(scRNA4SeJeAPmyoFBEC@meta.data), var = "rowname")
scRNA4SeJeAPmyoFBEC_tibble$rowname=rownames(scRNA4SeJeAPmyoFBEC@meta.data)
# 使用left_join合并数据框，然后根据结果替换df2中的celltype2列
merged_df <- left_join(scRNA4SeJeAPmyoFBEC_tibble, scRNA4SeJeAP_tibble, by = "rowname")
 
# 检查是否有NA值，如果有，则用原始的celltype2值填充
merged_df$cellType2 <- as.character(merged_df$cellType2)
merged_df$cellType.x<- as.character(merged_df$cellType.x)
# 使用ifelse进行替换
merged_df$celltype2 <- ifelse(is.na(merged_df$cellType2), merged_df$cellType.x, merged_df$cellType2)
 
# 如果需要，将celltype2转换回因子（并更新水平集）
# 注意：这一步只在celltype2原本就是因子且你需要保持它为因子类型时才需要
# levels_to_use <- unique(c(as.character(merged_df$celltype2), levels(as.factor(scRNASeJeAP@meta.data$cellType2))))
merged_df$celltype2 <- factor(merged_df$celltype2, levels = c('AP6','APNo6','myoFB','EC'))
scRNA4SeJeAPmyoFBEC@meta.data$celltype2=merged_df$celltype2

dot=DotPlot(scRNA4SeJeAPmyoFBEC, features = c(genes),group.by='celltype2')
dotdata=dot$data
write.csv('scRNA4SeJeAPmyoFBEC_expression.csv',dotdata)


# 关键基因散点图
library(readxl)
library(ggplot2)
library(ggrepel)

# 读取Excel文件
df <- read_excel("gene93.xlsx")
df <- df[,2:4]

# 绘制基本散点图（颜色按条件动态设置）
p <- ggplot(df, aes(x = `pct.exp-AP6`, y = perFCmyoFBvsAP6)) +
  # 核心修改：按条件设定颜色
  geom_point(
    aes(color = (`pct.exp-AP6` > 15 & perFCmyoFBvsAP6 > 3)), # 颜色映射逻辑
    size = 3, 
    alpha = 0.7
  ) +
  scale_color_manual(
    values = c("FALSE" = "steelblue", "TRUE" = "red") # 颜色映射规则
  ) +
  guides(color = "none") +  # 隐藏颜色图例（按需保留或删除此行）
  
  # 参考线
  geom_vline(xintercept = 15,lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = 3,lty=4,col="black",lwd=0.8) +
  
  # 坐标轴和标题
  labs(
    x = "pct.exp-AP6 (%)", 
    y = "perFC(myoFB vs AP6)",
    title = ""
  ) +
  scale_x_continuous(limits = c(0, 40))+
  scale_y_continuous(limits = c(0,6))+
  
  # 主题美化
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(size = 12),
    panel.grid.major = element_line(color = "grey90")
  ) +
  
  geom_text_repel(
    data = df,
    aes(label = label),
    box.padding = 0.5,
    segment.color = "grey50",
    max.overlaps = 20
  )
ggsave(filename = 'point93gene2.pdf',plot=p,width = 6,height = 6)


# beige 脂肪细胞关键基因分析

# 根据棕色脂肪细胞AP6差异基因对beige脂肪细胞做富集评分，找到类似AP6细胞后再进行关键基因的查找
# 提取差异基因（示例）
degs_bat <- FindMarkers(
  scRNA4SeJeAP,
  ident.1 = "6",  # 目标亚群
  ident.2 = NULL,           # 可设为其他亚群或留空
  min.pct = 0.25,
  logfc.threshold = 0.5
)
# 筛选显著基因（按p值或logFC）
sig_genes <- rownames(subset(degs_bat, p_val_adj < 0.05 & avg_log2FC > 0.5))

scRNA4SiJi=subset(x = scRNA4, subset = (orig.ident == "Si"|orig.ident == "Ji"))
scRNA4SiJiAP=subset(x = scRNA4SiJi, subset = (cellType == "AP"))
# 提取表达矩阵（log-normalized或SCT）
expr_matrix <- GetAssayData(scRNA4SiJiAP, slot = "data")  # log-normalized推荐

# 1. 构建基因集（列表形式）
gene_sets <- list(BAT_signature = sig_genes)  # 可添加多个基因集

# 2. 计算细胞级别AUC分数
cells_rankings <- AUCell_buildRankings(
  expr_matrix,
  nCores = 50,      # 多核加速
  plotStats = TRUE # 检查表达分布
)

cells_AUC <- AUCell_calcAUC(
  gene_sets,
  cells_rankings,
  aucMaxRank = ceiling(nrow(cells_rankings))  # 默认取前5%的基因计算AUC
)

# 3. 提取AUC分数并添加到Seurat对象
auc_scores <- getAUC(cells_AUC)
rownames(auc_scores) <- names(gene_sets) 
scRNA4SiJiAP[["BAT_AUC"]] <- CreateAssayObject(data = auc_scores)  # 作为新Assay存储

scRNA4SiJiAP@assays$BAT_AUC@key <- "rna_"
DefaultAssay(scRNA4SiJiAP) <- "BAT_AUC"
scRNA4SiJiAP@meta.data$orig.ident=factor(scRNA4SiJiAP@meta.data$orig.ident,levels =c('Ji','Si'))

# 将AUC分数添加到metadata（方便绘图）
scRNA4SiJiAP$BAT_AUC_score <- auc_scores["BAT_signature", ]

# 定义高AUC细胞（如Top 20%）
high_auc_cells <- colnames(scRNA4SiJiAP)[scRNA4SiJiAP$BAT_AUC_score > quantile(scRNA4SiJiAP$BAT_AUC_score, 0.8)]
scRNA4SiJiAP$BAT_like <- ifelse(colnames(scRNA4SiJiAP) %in% high_auc_cells, "BAT-like", "Other")

# AP亚群分类
DefaultAssay(scRNA4SiJiAP) <- "RNA"
# DefaultAssay(scRNA4SiJiAP) <- "integrated"

scRNA4SiJiAP <- NormalizeData(scRNA4SiJiAP, normalization.method = "LogNormalize", scale.factor = 10000)

scRNA4SiJiAP <- FindVariableFeatures(scRNA4SiJiAP,nfeatures = 3000)
scRNA4SiJiAP <- ScaleData(scRNA4SiJiAP,vars.to.regress = 'percent.mt')
scRNA4SiJiAP <- RunPCA(scRNA4SiJiAP)
scRNA4SiJiAP <- FindNeighbors(scRNA4SiJiAP, reduction = "pca", dims = 1:30)
scRNA4SiJiAP <- FindClusters(scRNA4SiJiAP, resolution = 0.3)
scRNA4SiJiAP <- RunUMAP(scRNA4SiJiAP, reduction = "pca", dims = 1:30) 
library(harmony)
scRNA4SiJiAP <- RunHarmony(scRNA4SiJiAP, group.by.vars = "orig.ident")
scRNA4SiJiAP  <- FindNeighbors(scRNA4SiJiAP , dims = 1:30 , reduction = "harmony")
scRNA4SiJiAP  <- FindClusters(scRNA4SiJiAP, save.snn=T , resolution = 0.3)
scRNA4SiJiAP  <- RunUMAP(scRNA4SiJiAP , dims=1:30,reduction='harmony')

saveRDS(scRNA4SiJiAP,'scRNA4SiJiAP_harmony.rds')

DefaultAssay(scRNA4SiJiAP) <- "BAT_AUC"
pdf(file=sprintf("scRNA4SiJiAP_genesets.AUCdotplot.pdf"),height = 3, width = 5)
DotPlot(scRNA4SiJiAP,group.by = 'seurat_clusters', features =c("BAT-signature")) + coord_flip()
dev.off()


# 获取不同类型细胞所有基因的表达平均值及表达比例
scRNA4SiJiAPmyoFBEC=subset(x = scRNA4SiJi, subset = (cellType == "AP"|cellType == "myoFB"|cellType == "EC"))

scRNA4SiJiAP@meta.data$celltype2<- 'APNo4'
Idents(scRNA4SiJiAP)='seurat_clusters'
filtered_cellsAP4 <- WhichCells(scRNA4SiJiAP,idents='4')
Idents(scRNA4SiJiAP)='celltype2'
Idents(scRNA4SiJiAP, cells = filtered_cellsAP4) <- 'AP4'
table(Idents(scRNA4SiJiAP))
scRNA4SiJiAP@meta.data$celltype2=Idents(scRNA4SiJiAP)

library(dplyr)
 
# 将df1转换为tibble，以便使用left_join，并设置行名为一个临时列
scRNA4SiJiAP_tibble <- tibble::rownames_to_column(as_tibble(scRNA4SiJiAP@meta.data), var = "rowname")
scRNA4SiJiAP_tibble$rowname=rownames(scRNA4SiJiAP@meta.data)
scRNA4SiJiAPmyoFBEC_tibble <- tibble::rownames_to_column(as_tibble(scRNA4SiJiAPmyoFBEC@meta.data), var = "rowname")
scRNA4SiJiAPmyoFBEC_tibble$rowname=rownames(scRNA4SiJiAPmyoFBEC@meta.data)
# 使用left_join合并数据框，然后根据结果替换df2中的celltype2列
merged_df <- left_join(scRNA4SiJiAPmyoFBEC_tibble, scRNA4SiJiAP_tibble, by = "rowname")
 
# 检查是否有NA值，如果有，则用原始的celltype2值填充
merged_df$celltype2 <- as.character(merged_df$celltype2)
merged_df$cellType.x<- as.character(merged_df$cellType.x)
# 使用ifelse进行替换
merged_df$celltype_f <- ifelse(is.na(merged_df$celltype2), merged_df$cellType.x, merged_df$celltype2)
 
# 如果需要，将celltype2转换回因子（并更新水平集）
# 注意：这一步只在celltype2原本就是因子且你需要保持它为因子类型时才需要
# levels_to_use <- unique(c(as.character(merged_df$celltype2), levels(as.factor(scRNASeJeAP@meta.data$cellType2))))
merged_df$celltype_f<- factor(merged_df$celltype_f, levels = c('AP4','APNo4','myoFB','EC'))
scRNA4SiJiAPmyoFBEC@meta.data$celltype2=merged_df$celltype_f
table(scRNA4SiJiAPmyoFBEC@meta.data$celltype2)
saveRDS(scRNA4SiJiAPmyoFBEC,'scRNA4SiJiAPmyoFBEC.rds')

library(limma)
library(edgeR)
data=scRNA4SiJiAP@assays$RNA@data

dd_CM<- as.data.frame(as.matrix(data))
# colnames(dd_CM)

group2<-  factor(scRNA4SiJiAP@meta.data$celltype2,levels = c("APNo4", "AP4"))#顺序靠后的Se vs 顺序靠前的Je;#对照组在前，处理组在后
table(group2)
dgelist <- DGEList(counts =dd_CM, group = group2)
#（2）过滤 low count 数据，例如 CPM 标准化（推荐）
keep <- rowSums(cpm(dgelist) > 1 ) >= 2
dgelist <- dgelist[keep, , keep.lib.sizes = FALSE]
#（3）标准化，以 TMM 标准化为例
dgelist_norm <- calcNormFactors(dgelist, method = 'TMM')
#差异表达基因分析
design <- model.matrix(~group2)
#（1）估算基因表达值的离散度
dge <- estimateDisp(dgelist_norm, design, robust = TRUE)
#（2）模型拟合，edgeR 提供了多种拟合算法
#负二项广义对数线性模型
fit <- glmFit(dge, design, robust = TRUE)
lrt <- topTags(glmLRT(fit), n = nrow(dgelist$counts))
sample1='SiJiAP4'
sample2='SiJiAPNo4'
write.table(lrt, sprintf('%svs%s.txt',sample1,sample2), sep = '\t', col.names = NA, quote = FALSE)
library(ggplot2)

d1<- read.table(sprintf('%svs%s.txt',sample1,sample2))
colnames(d1)
d1$gene<- rownames(d1)

# 设置p_value和logFC的阈值
cut_off_FDR = 0.05  #统计显著性
cut_off_logFC = 0.5           #差异倍数值
d1$change<- 'Stable'
# 根据阈值参数，上调基因设置为‘up’，下调基因设置为‘Down’，无差异设置为‘Stable’，并保存到change列中
d1$change = ifelse(d1$FDR< cut_off_FDR & abs(d1$logFC) > cut_off_logFC, 
                        ifelse(d1$logFC> cut_off_logFC ,'Up','Down'),
                        'Stable')
head(d1)
table(d1$change)
write.csv(file='sijiAP4diff.csv',d1)
p <- ggplot(
  # 数据、映射、颜色
  d1, aes(x = logFC, y = -log10(FDR), colour=change)) +
  geom_point(alpha=0.4, size=2) +
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
  # 辅助线
  geom_vline(xintercept=c(-0.5,0.5),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_FDR),lty=4,col="black",lwd=0.8) +
  # 坐标轴
  labs(x="log2(fold change)",
       y="-log10 (FDR)")+
  scale_x_continuous(limits = c(-2.5, 2.5))+
  scale_y_continuous(limits = c(0,10))+
  theme_bw()+
  # 图例
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())+
  ggtitle('AP4 v.s. APNo4')
# p
ggsave(filename = "sijidiffgene2nolabelAP4.pdf", height = 6, width = 6, plot = p)

genes=rownames(d1)[d1$change=='Up']

DefaultAssay(scRNA4SiJiAPmyoFBEC)='RNA'
dot=DotPlot(scRNA4SiJiAPmyoFBEC, features = c(genes),group.by='celltype2')
dotdata=dot$data
write.csv('scRNA4SiJiAPmyoFBEC_expression.csv',dotdata)



library(readxl)
library(ggplot2)
library(ggrepel)

# 读取Excel文件
data <- read_excel("gene49.xlsx")
data <- data[,2:4]

df <- data
# 绘制基础散点图


# 绘制基本散点图（颜色按条件动态设置）
p <- ggplot(df, aes(x = `pct.exp-AP4`, y = perFCmyoFBvsAP4)) +
  # 核心修改：按条件设定颜色
  geom_point(
    aes(color = (`pct.exp-AP4` > 15 & perFCmyoFBvsAP4 > 3)), # 颜色映射逻辑
    size = 3, 
    alpha = 0.7
  ) +
  scale_color_manual(
    values = c("FALSE" = "steelblue", "TRUE" = "red") # 颜色映射规则
  ) +
  guides(color = "none") +  # 隐藏颜色图例（按需保留或删除此行）
  
  # 参考线
  geom_vline(xintercept = 15,lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = 3,lty=4,col="black",lwd=0.8) +
  
  # 坐标轴和标题
  labs(
    x = "pct.exp-AP4 (%)", 
    y = "perFC(myoFB vs AP4)",
    title = ""
  ) +
  scale_x_continuous(limits = c(0, 40))+
  scale_y_continuous(limits = c(0,6))+
  
  # 主题美化
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(size = 12),
    panel.grid.major = element_line(color = "grey90")
  ) +
  
  geom_text_repel(
    data = df,
    aes(label = label),
    box.padding = 0.5,
    segment.color = "grey50",
    max.overlaps = 20
  )
ggsave(filename = 'point49gene.pdf',plot=p,width = 6,height = 6)