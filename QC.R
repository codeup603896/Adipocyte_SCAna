library(Seurat)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(dplyr)
library(tidyr)


# 样本matrix数据路径
dir4 = c(
  "../Je_Ana/step3/filtered_feature_bc_matrix",
  "../Ji_Ana/step3/filtered_feature_bc_matrix",
  "../Se_Ana/step3/filtered_feature_bc_matrix",
  "../Si_Ana/step3/filtered_feature_bc_matrix"
)

scRNAlist4 <- list()  # 创建一个空的列表
names(dir4) = c("Je", "Ji",'Se','Si')
sample_name = c("Je", "Ji",'Se','Si')

for (i in 1:length(dir4)){
  counts <- Read10X(data.dir = dir4[i])

  scRNAlist4[[i]] <- CreateSeuratObject(counts = counts,min.cells = 3,min.features = 200)
  scRNAlist4[[i]] <- RenameCells(scRNAlist4[[i]],add.cell.id = sample_name[i])
}

# 合并Seurat object
seurat_obj <- merge(scRNAlist4[[1]],y = c(scRNAlist4[[2]],scRNAlist4[[3]],scRNAlist4[[4]]))

# 添加线粒体基因百分比计算 (假设是人数据)
# seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# 如果是小鼠数据，使用:
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")

# 设置QC阈值
min_genes <- 200     # 每个核至少检测到的基因数
max_genes <- 4000    # 每个核最多检测到的基因数(过滤可能的双联体)
max_mt <- 5        # 最大线粒体基因百分比

# 查看原始QC指标
head(seurat_obj@meta.data)

# 绘制QC指标分布图
qc_plots <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                    ncol = 3, group.by = "orig.ident") & 
  theme(plot.title = element_text(size=10))
# qc_plots
ggsave(filename = "qc_plots.pdf", height = 6, width = 6, plot = qc_plots)
# 检查基因数与UMI数的关系
plot1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  geom_hline(yintercept = max_mt, linetype = "dashed", color = "red")
ggsave(filename = "Count_mt_cor.pdf", height = 6, width = 6, plot =plot1)
plot2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_hline(yintercept = min_genes, linetype = "dashed", color = "red") +
  geom_hline(yintercept = max_genes, linetype = "dashed", color = "red")
ggsave(filename = "Count_Feature_cor.pdf", height = 6, width = 6, plot =plot2)

# 显示过滤前每个样本的细胞计数
pre_filter_counts <- table(seurat_obj$orig.ident)
cat("过滤前的细胞计数:\n")
print(pre_filter_counts)

# 实施QC过滤
seurat_obj_filtered <- subset(seurat_obj, 
                              subset = nFeature_RNA > min_genes & 
                                nFeature_RNA < max_genes & 
                                percent.mt < max_mt)

# 显示过滤后每个样本的细胞计数
post_filter_counts <- table(seurat_obj_filtered$orig.ident)
cat("\n过滤后的细胞计数:\n")
print(post_filter_counts)

# 计算过滤百分比
filter_stats <- data.frame(
  Sample = names(pre_filter_counts),
  Before_Filter = as.numeric(pre_filter_counts),
  After_Filter = as.numeric(post_filter_counts[match(names(pre_filter_counts), names(post_filter_counts))])
)

filter_stats$Percentage_Remaining <- round(filter_stats$After_Filter / filter_stats$Before_Filter * 100, 2)
filter_stats$Cells_Removed <- filter_stats$Before_Filter - filter_stats$After_Filter
filter_stats$Percentage_Removed <- 100 - filter_stats$Percentage_Remaining

cat("\n过滤统计:\n")
print(filter_stats)

# 绘制过滤前后细胞计数对比图

filter_stats_long <- gather(filter_stats, key = "Filter_Status", value = "Cell_Count", Before_Filter, After_Filter)

p_compare=ggplot(filter_stats_long, aes(x = Sample, y = Cell_Count, fill = Filter_Status)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(aes(label = Cell_Count), vjust = -0.3, position = position_dodge(0.9), size = 3) +
  labs(title = "Comparison", y = "cell num", x = "samples") +
  theme_bw() +theme(panel.border = element_blank(), axis.line = element_line(),panel.grid=element_blank()) +
scale_y_continuous(expand = c(0,0),limits = c(0,12000)) 
ggsave(filename = "p_compare.pdf", height = 4, width =6, plot =p_compare)

# 绘制过滤后QC指标
p_qc_filter=VlnPlot(seurat_obj_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, group.by = "orig.ident") & 
  theme(plot.title = element_text(size=10))
ggsave(filename = "p_qc_filter.pdf", height = 6, width = 6, plot =p_qc_filter)
